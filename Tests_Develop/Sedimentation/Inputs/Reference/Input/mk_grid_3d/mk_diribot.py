import os
from glob import glob
import sys

# === 1. Leggi nbot da grid3D ===
with open('grid3D', 'r') as f:
    f.readline()  # salta la prima riga
    seconda_riga = f.readline()
    nbot = int(seconda_riga.strip().split()[0])

# === 2. Apri file di output ===
with open('diri_bot', 'w') as f_out:
    # Scrivi nbot
    f_out.write(f"{nbot}\n")

    # Scrivi numeri da 1 a nbot, 6 per riga
    for i in range(1, nbot + 1):
        f_out.write(f"{i:15d} ")
        if i % 6 == 0 or i == nbot:
            f_out.write("\n")

    # === 3. Leggi tutti i file di dati mensili ===
    folder = 'hydro_step01'
    pattern = os.path.join(folder, 'MULTIPLY_HEAD_*.txt')
    files = sorted(glob(pattern))  # ordina per data (importante!)
    # Controllo: se non trova file, esce con messaggio
    if not files:
       print(f"Errore: nessun file trovato nel percorso '{folder}' con pattern '{pattern}'")
       sys.exit(1)
    

    # Tempo iniziale
    dt = 1.0 / 12  # un mese in frazioni di anno
    time = 100.0
    
    f_out.write(f"{time:.3f}     !time0\n")
    for file_path in files:
        f_out.write(f"{time:.3f}     !time\n")
        with open(file_path, 'r') as f_in:
            next(f_in)  # salta intestazione

            valori = []
            for line in f_in:
                parts = line.strip().split()
                if len(parts) >= 3:
                    valori.append(parts[2])  # terzo valore

            # Scrivi valori, 6 per riga
            for i, val in enumerate(valori, 1):
                f_out.write(f"{float(val):14.6e} ")
                if i % 6 == 0 or i == len(valori):
                    f_out.write("\n")

        time += dt  # incrementa il tempo per il file successivo
