# Function to read particle distribution from file entered by user

def get_file_selection(h5_files, files_dir):
    # Display the list of .h5 files
    print("Available simulation .h5 files:")
    for i, file in enumerate(h5_files):
        print(f"{i+1}. {file}")

    # Prompt the user to select a file
    selection = input("Enter the number corresponding to the file you want to select: ")

    # Validate the user's input and retrieve the selected filepath
    try:
        selection = int(selection)
        if 1 <= selection <= len(h5_files):
            selected_file = os.path.join(files_dir, h5_files[selection-1])
            print(f"You selected: {selected_file}")
            return selected_file
        else:
            print("Invalid selection.")
    except ValueError:
        print("Invalid input. Please enter a number.")
        
