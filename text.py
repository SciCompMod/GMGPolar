def read_and_format_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Find the indices where the data starts
    original_index = lines.index('Original Extrapolated Smoother\n') + 1
    refactored_index = lines.index('Refactored Extrapolated Smoother\n') + 1
    
    # Extract data lines
    original_data = lines[original_index].strip().split()
    refactored_data = lines[refactored_index].strip().split()
    
    # Prepare the header
    header = "Original Extrapolated Smoother  Refactored Extrapolated Smoother"
    
    # Combine the data into a formatted string
    formatted_lines = [header]
    for original, refactored in zip(original_data, refactored_data):
        formatted_lines.append(f"{original}  {refactored}")
    
    # Print the formatted table
    for line in formatted_lines:
        print(line)

    # Optionally, save the table to a file
    with open('output_table.txt', 'w') as output_file:
        for line in formatted_lines:
            output_file.write(line + '\n')

# Replace 'input.txt' with the path to your actual text file
read_and_format_file('text.txt')
