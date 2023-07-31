import tkinter as tk

def calculate_collimator_focus():
    try:
        temp_mount = float(temp_mount_entry.get())
        tilt = float(tilt_entry.get())
        obf = obf_var.get()

        if obf == 'Yes':
            obf = 1
        else:
            obf = 0

        cf = 11.0 - 0.08 * temp_mount - 0.14 * (tilt - 25) + 0.7 * obf

        result_label.config(text=f'Collimator focus calculated at {cf}mm')

    except ValueError:
        result_label.config(text='Invalid input. Please enter numeric values for temperature and tilt.')

root = tk.Tk()
root.title('DeVeny Collimator Focus Calculator')

# Create input fields
temp_mount_label = tk.Label(root, text='Mount temperature in deg C:')
temp_mount_label.pack()
temp_mount_entry = tk.Entry(root)
temp_mount_entry.pack()

tilt_label = tk.Label(root, text='Grating tilt in degrees:')
tilt_label.pack()
tilt_entry = tk.Entry(root)
tilt_entry.pack()

obf_label = tk.Label(root, text='Is there an order blocking filter present?')
obf_label.pack()

obf_var = tk.StringVar(root, value='No')
obf_yes_radio = tk.Radiobutton(root, text='Yes', variable=obf_var, value='Yes')
obf_yes_radio.pack()
obf_no_radio = tk.Radiobutton(root, text='No', variable=obf_var, value='No')
obf_no_radio.pack()

calculate_button = tk.Button(root, text='Calculate', command=calculate_collimator_focus)
calculate_button.pack()

result_label = tk.Label(root, text='')
result_label.pack()

root.mainloop()
