import os
import threading
import subprocess
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

def run_pipeline(input_file, ids_file, output_dir, reject, single_fasta, progress_bar):

    # Start the progress bar
    progress_bar.start()

    # Change working directory
    os.chdir(output_dir)
    
    # Convert windows to wsl pathways for the ids file
    ids_file_fixed = str(ids_file).replace(" ","\ ")

    # Convert windows to wsl pathways for the input fasta file
    fasta_file_fixed = str(input_file).replace(" ","\ ")

    # Output file
    if not single_fasta:
        if not reject:
            output_file = f"{os.path.splitext(os.path.basename(ids_file))[0]}.fasta"
        else:
            output_file = f"non_{os.path.splitext(os.path.basename(ids_file))[0]}.fasta"

        output_file_fixed = str(output_file).replace(" ","\ ")
    
    # Run command
    if not single_fasta:
        if not reject:
            command = f"fastaselecth -com -in {fasta_file_fixed} -sel {ids_file_fixed} -out {output_file_fixed}"
        else:
            command = f"fastaselecth -com -reject -in {fasta_file_fixed} -sel {ids_file_fixed} -out {output_file_fixed}"
    else:
            command = f"paste {ids_file_fixed} {ids_file_fixed} | fastaselecth -com -fragc -in {fasta_file_fixed} -sel \"-\" -out \"%s.fasta\""

    try:
        subprocess.run(["bash", "-c", command], check=True)
        progress_bar.stop()
        if not single_fasta:
            messagebox.showinfo("Success", f"Output file created at {os.path.abspath(output_file)}")
        else:
            messagebox.showinfo("Success", f"Output files created at {output_dir}")

    except subprocess.CalledProcessError as e:
        progress_bar.stop()
        messagebox.showerror(f"Error: {e}")
        
def start_thread():
    input_file = input_file_var.get()
    ids_file = ids_file_var.get()
    output_dir = output_dir_var.get()
    reject = reject_var.get()
    single_fasta = single_fasta_var.get()
    
    if not input_file:
        messagebox.showwarning("Input Error", "Please select an input FASTA file.")
        return

    if not ids_file:
        messagebox.showwarning("Input Error", "Please select an input TXT file.")
        return
    
    if not  output_dir:
        messagebox.showwarning("Input Error", "Please select an output directory.")
        return
    if single_fasta and reject:
        messagebox.showwarning("Input Error", "The reject option cannot be used with single-fasta files as output.")
        return        

    # Start command in a new thread
    thread = threading.Thread(target=run_pipeline, args=(input_file, ids_file, output_dir, reject, single_fasta, progress_bar))
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

# Checkbox for additional option
reject_var = tk.BooleanVar(value=False)
tk.Checkbutton(app, text=" Reject selected entries", variable=reject_var).grid(row=3, column=1, padx=10, pady=10, sticky="w")

# Checkbox for additional option
single_fasta_var = tk.BooleanVar(value=False)
tk.Checkbutton(app, text="Export to single-fasta files", variable=single_fasta_var).grid(row=4, column=1, padx=10, pady=10, sticky="w")

# Progress Bar (indeterminate)
progress_bar = ttk.Progressbar(app, mode="indeterminate", length=200)
progress_bar.grid(row=5, column=0, columnspan=3, padx=10, pady=20)

# Start button
tk.Button(app, text="Run program", command=start_thread).grid(row=6, column=1, padx=10, pady=20)

app.mainloop()
