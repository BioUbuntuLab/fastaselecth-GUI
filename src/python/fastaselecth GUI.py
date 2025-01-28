import os
import threading
import subprocess
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

def run_pipeline(input_file, ids_file, output_dir, progress_bar):

    # Start the progress bar
    progress_bar.start()

    # Change working directory
    os.chdir(output_dir)
    
    # Convert windows to wsl pathways for the ids file
    splitted_ids_pathway = str(ids_file).split(":",1)
    drive_ids = str(splitted_ids_pathway[0]).lower()
    ids_file_fixed = f"/mnt/{drive_ids}{splitted_ids_pathway[1]}".replace(" ","\ ")

    # Convert windows to wsl pathways for the input fasta file
    splitted_fasta_pathway = str(input_file).split(":",1)
    drive_fasta = str(splitted_fasta_pathway[0]).lower()
    fasta_file_fixed = f"/mnt/{drive_fasta}{splitted_fasta_pathway[1]}".replace(" ","\ ")

    # Output folder
    output_file = f"{os.path.splitext(os.path.basename(ids_file))[0]}.fasta"
    output_file_fixed = str(output_file).replace(" ","\ ")
    
    # Run command
    command = f"fastaselecth -in {fasta_file_fixed} -sel {ids_file_fixed} -out {output_file_fixed}"

    try:
        subprocess.run(["wsl", "bash", "-c", command], check=True, creationflags=subprocess.CREATE_NO_WINDOW)
        progress_bar.stop()
        messagebox.showinfo("Success", f"Output file created at {os.path.abspath(output_file)}")


    except subprocess.CalledProcessError as e:
        progress_bar.stop()
        messagebox.showerror(f"Error: {e}")
        
def start_thread():
    input_file = input_file_var.get()
    ids_file = ids_file_var.get()
    output_dir = output_dir_var.get()

    if not input_file:
        messagebox.showwarning("Input Error", "Please select an input FASTA file.")
        return

    if not ids_file:
        messagebox.showwarning("Input Error", "Please select an input TXT file.")
        return
    
    if not  output_dir:
        messagebox.showwarning("Input Error", "Please select an output directory.")
        return
    

    # Start command in a new thread
    thread = threading.Thread(target=run_pipeline, args=(input_file, ids_file, output_dir, progress_bar))
    thread.start()

def select_fasta_file():
    file_path = filedialog.askopenfilename()
    input_file_var.set(file_path)

def select_ids_file():
    file_path = filedialog.askopenfilename()
    ids_file_var.set(file_path)

def select_directory():
    file_path = filedialog.askdirectory()
    output_dir_var.set(file_path)
# Set up tkinter app
app = tk.Tk()
app.title("fastaselecth GUI")

# Input file selection
input_file_var = tk.StringVar()
tk.Label(app, text="Input FASTA File:").grid(row=0, column=0, padx=10, pady=10, sticky="e")
tk.Entry(app, textvariable=input_file_var, width=40).grid(row=0, column=1, padx=10, pady=10)
tk.Button(app, text="Browse", command=select_fasta_file).grid(row=0, column=2, padx=10, pady=10)

# Input file selection
ids_file_var = tk.StringVar()
tk.Label(app, text="Input 1-column TXT File With IDs:").grid(row=1, column=0, padx=10, pady=10, sticky="e")
tk.Entry(app, textvariable=ids_file_var, width=40).grid(row=1, column=1, padx=10, pady=10)
tk.Button(app, text="Browse", command=select_ids_file).grid(row=1, column=2, padx=10, pady=10)

# Input file selection
output_dir_var = tk.StringVar()
tk.Label(app, text="Output Directory:").grid(row=2, column=0, padx=10, pady=10, sticky="e")
tk.Entry(app, textvariable=output_dir_var, width=40).grid(row=2, column=1, padx=10, pady=10)
tk.Button(app, text="Browse", command=select_directory).grid(row=2, column=2, padx=10, pady=10)

# Progress Bar (indeterminate)
progress_bar = ttk.Progressbar(app, mode="indeterminate", length=200)
progress_bar.grid(row=3, column=0, columnspan=3, padx=10, pady=20)

# Start button
tk.Button(app, text="Run program", command=start_thread).grid(row=4, column=1, padx=10, pady=20)

app.mainloop()
