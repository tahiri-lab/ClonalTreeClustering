"""Interactive TUI for running the WMFD pipeline scripts (2 → 7.c) on
user-provided weighted Newick trees.

Left column  → tree list + multiline input + Add/Remove/Clear.
Right column → script list + Run.
After a script finishes, a modal shows its output file (truncated)
with "Volver" and "Ejecutar siguiente" actions.

Run with:  uv run python src/tui.py   (from the WMFD/ directory)
"""
from __future__ import annotations

import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from textual import work
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical
from textual.screen import ModalScreen
from textual.widgets import (
    Button,
    Footer,
    Header,
    Label,
    ListItem,
    ListView,
    Static,
    TextArea,
)


WMFD_DIR = Path(__file__).resolve().parent.parent
SRC_DIR = WMFD_DIR / "src"
DATA_DIR = WMFD_DIR / "simulated_data"
INPUT_FILE = DATA_DIR / ".tui_input.txt"
MAX_LINES = 500


@dataclass(frozen=True)
class ScriptEntry:
    filename: str
    output: str


SCRIPTS: list[ScriptEntry] = [
    ScriptEntry("2.Penalty.py", "penalties.txt"),
    ScriptEntry(
        "3.a.weight_BL_matrices_with_normalized_matrices_text_output.py",
        "matrices_with_normalization_60.txt",
    ),
    ScriptEntry("3.b.BL_W_sum_common_uncommon.py", "node_comparison_results.txt"),
    ScriptEntry(
        "3.c.BL_W_normalized_sum_common_uncommon.py",
        "node_comparison_normalized_results.txt",
    ),
    ScriptEntry("4.BL_W_differences_matrix.py", "normalized_pairwise_differences_60.txt"),
    ScriptEntry("5.b.adjacency_matrices.py", "adjacency_matrices.txt"),
    ScriptEntry("5.c.normalized_HD.py", "Normalized_HD.txt"),
    ScriptEntry("6.a.height.py", "height_matrices_60.txt"),
    ScriptEntry("6.b.height_sum_commom_uncommon.py", "height_sums_common_uncommon.txt"),
    ScriptEntry(
        "6.c.height_normalized_sum_common_uncommon.py",
        "height_sums_common_uncommon_normalized.txt",
    ),
    ScriptEntry("7.a.degree.py", "tree_degrees_results.txt"),
    ScriptEntry("7.b.degree_sum_common_uncommon.py", "degree_comparison_results.txt"),
    ScriptEntry(
        "7.c.degree_normalized_sum_common_uncommon.py",
        "degree_comparison_normalized_results.txt",
    ),
]

MAX_SEQ_NUM_HARDCODED = 28
HARDCODED_SCRIPTS = ("5.b.adjacency_matrices.py", "5.c.normalized_HD.py")

_PREFIX_RE = re.compile(r"^[^(]*?:\s")


def _normalize_tree(line: str, idx: int) -> str:
    line = line.strip()
    if not line:
        return ""
    if _PREFIX_RE.match(line):
        return line
    return f"tui_{idx + 1}: {line}"


def _max_seq_num(trees: list[str]) -> int:
    nums = [int(m) for tree in trees for m in re.findall(r"seq(\d+)", tree)]
    return max(nums) if nums else 0


def _truncate(text: str, max_lines: int) -> tuple[str, int]:
    lines = text.splitlines()
    if len(lines) <= max_lines:
        return text, 0
    head = lines[:max_lines]
    overflow = len(lines) - max_lines
    head.append(f"... [truncado, {overflow} líneas más]")
    return "\n".join(head), overflow


@dataclass
class RunResult:
    script: ScriptEntry
    returncode: int
    stdout: str
    stderr: str
    output_path: Path
    output_text: str
    output_size: int
    truncated_lines: int


class ResultScreen(ModalScreen[str]):
    """Shows the output of a script run."""

    BINDINGS = [
        Binding("escape", "back", "Volver"),
        Binding("n", "next", "Siguiente"),
    ]

    DEFAULT_CSS = """
    ResultScreen {
        align: center middle;
    }
    ResultScreen > Vertical {
        background: $surface;
        border: thick $primary;
        padding: 1 2;
        width: 95%;
        height: 90%;
    }
    ResultScreen #header {
        height: auto;
        padding: 0 0 1 0;
    }
    ResultScreen #body {
        height: 1fr;
    }
    ResultScreen #footer {
        height: 3;
        align: center middle;
    }
    ResultScreen #footer Button {
        margin: 0 1;
    }
    """

    def __init__(self, result: RunResult, has_next: bool) -> None:
        super().__init__()
        self.result = result
        self.has_next = has_next

    def compose(self) -> ComposeResult:
        r = self.result
        if r.returncode == 0:
            trunc = (
                f" — truncado ({r.truncated_lines} líneas más)"
                if r.truncated_lines
                else ""
            )
            header = (
                f"[b]{r.script.filename}[/]  →  exit 0\n"
                f"Salida: {r.output_path}  ({r.output_size} bytes){trunc}"
            )
            body_text = r.output_text or "(archivo vacío o inexistente)"
        else:
            header = (
                f"[b]{r.script.filename}[/]  →  exit {r.returncode}  [red]ERROR[/]\n"
                f"stderr:"
            )
            body_text = r.stderr or r.stdout or "(sin salida)"

        with Vertical():
            yield Static(header, id="header", markup=True)
            yield TextArea(body_text, id="body", read_only=True, soft_wrap=False)
            with Horizontal(id="footer"):
                yield Button("← Volver", id="back")
                yield Button(
                    "Ejecutar siguiente →",
                    id="next",
                    variant="primary",
                    disabled=(not self.has_next) or r.returncode != 0,
                )

    def on_button_pressed(self, event: Button.Pressed) -> None:
        self.dismiss(event.button.id)

    def action_back(self) -> None:
        self.dismiss("back")

    def action_next(self) -> None:
        if self.has_next and self.result.returncode == 0:
            self.dismiss("next")


class WMFDApp(App[None]):
    CSS = """
    #main {
        layout: horizontal;
        height: 1fr;
    }
    #left, #right {
        width: 1fr;
        border: solid $primary;
        padding: 1;
    }
    #left {
        margin-right: 1;
    }
    #trees, #scripts {
        height: 1fr;
        border: round $secondary;
    }
    #new_tree {
        height: 6;
        margin-top: 1;
    }
    #left_buttons, #right_buttons {
        height: 3;
        margin-top: 1;
        align: center middle;
    }
    #left_buttons Button, #right_buttons Button {
        margin: 0 1;
    }
    #status {
        height: 2;
        color: $warning;
        margin-top: 1;
    }
    """

    BINDINGS = [
        Binding("ctrl+q", "quit", "Salir"),
    ]

    def __init__(self) -> None:
        super().__init__()
        self.trees: list[str] = []
        self.last_script_idx: int | None = None

    def compose(self) -> ComposeResult:
        yield Header(show_clock=False)
        with Horizontal(id="main"):
            with Vertical(id="left"):
                yield Label("[b]Árboles actuales[/]", markup=True)
                yield ListView(id="trees")
                yield TextArea(
                    "",
                    id="new_tree",
                    soft_wrap=True,
                )
                yield Static("", id="status", markup=True)
                with Horizontal(id="left_buttons"):
                    yield Button("Agregar", id="add", variant="success")
                    yield Button("Eliminar", id="remove", variant="warning")
                    yield Button("Limpiar", id="clear", variant="error")
            with Vertical(id="right"):
                yield Label("[b]Scripts[/]", markup=True)
                yield ListView(
                    *(ListItem(Label(s.filename)) for s in SCRIPTS),
                    id="scripts",
                )
                with Horizontal(id="right_buttons"):
                    yield Button("Ejecutar seleccionado", id="run", variant="primary")
        yield Footer()

    def on_mount(self) -> None:
        self.title = "WMFD TUI"
        self.sub_title = "scripts 2 → 7.c"
        self.query_one("#scripts", ListView).index = 0
        self._refresh_status()

    # --- tree management ---

    def _refresh_trees_view(self) -> None:
        lv = self.query_one("#trees", ListView)
        lv.clear()
        for i, t in enumerate(self.trees):
            preview = (t[:80] + "…") if len(t) > 80 else t
            lv.append(ListItem(Label(f"{i + 1}. {preview}")))
        self._refresh_status()

    def _refresh_status(self) -> None:
        status = self.query_one("#status", Static)
        run_btn = self.query_one("#run", Button)
        msgs: list[str] = []
        if len(self.trees) < 2:
            run_btn.disabled = True
            msgs.append(f"Pegue al menos 2 árboles ({len(self.trees)}/2).")
        else:
            run_btn.disabled = False
            msgs.append(f"{len(self.trees)} árboles cargados.")
        ms = _max_seq_num(self.trees)
        if ms > MAX_SEQ_NUM_HARDCODED:
            msgs.append(
                f"[orange1]Aviso:[/] hay seq{ms}; "
                f"{', '.join(HARDCODED_SCRIPTS)} truncan en seq{MAX_SEQ_NUM_HARDCODED}."
            )
        status.update("  ".join(msgs))

    def action_add_tree(self) -> None:
        ta = self.query_one("#new_tree", TextArea)
        raw = ta.text.strip()
        if not raw:
            return
        for line in raw.splitlines():
            normalized = _normalize_tree(line, len(self.trees))
            if normalized:
                self.trees.append(normalized)
        ta.text = ""
        self._refresh_trees_view()

    def action_remove_tree(self) -> None:
        lv = self.query_one("#trees", ListView)
        idx = lv.index
        if idx is None or idx < 0 or idx >= len(self.trees):
            return
        del self.trees[idx]
        self._refresh_trees_view()

    def action_clear_trees(self) -> None:
        self.trees.clear()
        self._refresh_trees_view()

    # --- script execution ---

    def _selected_script_idx(self) -> int:
        lv = self.query_one("#scripts", ListView)
        return lv.index if lv.index is not None else 0

    def _run_at(self, idx: int) -> None:
        if idx < 0 or idx >= len(SCRIPTS):
            return
        if len(self.trees) < 2:
            self.bell()
            return
        self.last_script_idx = idx
        DATA_DIR.mkdir(exist_ok=True)
        INPUT_FILE.write_text("\n".join(self.trees) + "\n", encoding="utf-8")
        run_btn = self.query_one("#run", Button)
        run_btn.disabled = True
        run_btn.label = "Ejecutando..."
        self._spawn(idx)

    @work(thread=True, exclusive=True)
    def _spawn(self, idx: int) -> None:
        entry = SCRIPTS[idx]
        script_path = SRC_DIR / entry.filename
        try:
            proc = subprocess.run(
                [sys.executable, str(script_path), "--input", str(INPUT_FILE)],
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                cwd=str(SRC_DIR),
            )
            output_path = DATA_DIR / entry.output
            if proc.returncode == 0 and output_path.exists():
                raw = output_path.read_text(encoding="utf-8", errors="replace")
                txt, truncated = _truncate(raw, MAX_LINES)
                output_size = output_path.stat().st_size
            else:
                txt, truncated, output_size = "", 0, 0
            result = RunResult(
                script=entry,
                returncode=proc.returncode,
                stdout=proc.stdout,
                stderr=proc.stderr,
                output_path=output_path,
                output_text=txt,
                output_size=output_size,
                truncated_lines=truncated,
            )
        except Exception as exc:
            result = RunResult(
                script=entry,
                returncode=-1,
                stdout="",
                stderr=f"Excepción al ejecutar: {exc}",
                output_path=DATA_DIR / entry.output,
                output_text="",
                output_size=0,
                truncated_lines=0,
            )
        self.call_from_thread(self._show_result, result)

    def _show_result(self, result: RunResult) -> None:
        run_btn = self.query_one("#run", Button)
        run_btn.label = "Ejecutar seleccionado"
        self._refresh_status()
        idx = SCRIPTS.index(result.script)
        has_next = idx + 1 < len(SCRIPTS)
        self.push_screen(ResultScreen(result, has_next), self._on_result_dismiss)

    def _on_result_dismiss(self, action: str | None) -> None:
        if action == "next" and self.last_script_idx is not None:
            next_idx = self.last_script_idx + 1
            if next_idx < len(SCRIPTS):
                self.query_one("#scripts", ListView).index = next_idx
                self._run_at(next_idx)

    # --- buttons ---

    def on_button_pressed(self, event: Button.Pressed) -> None:
        bid = event.button.id
        if bid == "add":
            self.action_add_tree()
        elif bid == "remove":
            self.action_remove_tree()
        elif bid == "clear":
            self.action_clear_trees()
        elif bid == "run":
            self._run_at(self._selected_script_idx())


if __name__ == "__main__":
    WMFDApp().run()
