filename = input("Entrez le nom du fichier : ")
counts = [0] * 61

with open(filename, 'r') as file:
    for line in file:
        for i in range(1, 61):
            if ">seq" + str(i) in line:
                counts[i] += 1

total_count = sum(counts)
print("Le nombre total d'occurrences de '>seqN' dans le fichier est :", total_count)
for i in range(1, 61):
    print("Le nombre d'occurrences de '>seq" + str(i) + "' dans le fichier est :", counts[i])