def read_newick(file_path):
    """
    Reads a Newick-formatted string from a file, ensuring it ends with a semicolon.

    Parameters:
    -----------
    file_path : str
        Path to the file containing the Newick string.

    Returns:
    --------
    newick : str
        The Newick string read from the file.

    Raises:
    -------
    ValueError
        If the end of file is reached before encountering a terminating semicolon (';').

    Side effects:
    -------------
    Runs a format check using NewickCheck().
    """
    with open(file_path, 'r') as f:
        newick = ''
        while True:
            c = f.read(1)
            if not c:
                raise ValueError("Unexpected end of file before finding ';'")
            newick += c
            if c == ';':
                break

    NewickCheck(newick)  # Validate format
    return newick


def NewickCheck(Newick):
    """
    Verifies the basic syntactic validity of a Newick string.

    Parameters:
    -----------
    Newick : str
        The Newick-formatted tree string to be validated.

    Behavior:
    ---------
    - Prints error messages for various common formatting mistakes.
    - Returns None after each message.
    - Prints "correct format of Newick" if all checks pass.

    Checks performed:
    -----------------
    - Starts with an opening parenthesis '('.
    - Balanced number of opening and closing parentheses.
    - Contains exactly one terminating semicolon ';'.
    - Disallows '%' character, which is not valid in Newick.
    - (Optional edge length check via ':' is commented out)
    """

    colon_count = 0
    paren_balance = 0
    semicolon_count = 0
    percent_count = 0

    for char in Newick:
        if char == ':':
            colon_count += 1
        elif char == '(':
            paren_balance += 1
        elif char == ')':
            paren_balance -= 1
        elif char == ';':
            semicolon_count += 1
        elif char == '%':
            percent_count += 1

    if Newick[0] != '(':
        print("Incorrect Newick file format. Newick string must begin with a '(' character.")
        return

    # Optional: Uncomment to require edge lengths
    # if colon_count == 0:
    #     print("Incorrect Newick file format. Edge lengths must be indicated after ':' characters.")
    #     return

    if paren_balance < 0:
        print("Incorrect Newick file format. Missing opening '(' character(s).")
        return
    elif paren_balance > 0:
        print("Incorrect Newick file format. Missing closing ')' character(s).")
        return

    if semicolon_count == 0:
        print("Incorrect Newick file format. Newick string must be terminated with a single ';' character.")
        return
    elif semicolon_count > 1:
        print("Incorrect Newick file format. Only one semicolon ';' should appear at the end.")
        return

    if percent_count > 0:
        print("Incorrect Newick file format. '%' character is not allowed.")
        return

    print("Correct format of Newick.")
    return
