import re

# Define the filename
filename = 'massif.out.2908275'

# Initialize variables to store maximum memory usage
max_memory_usage_MB = 0

# Use regular expressions to find all occurrences of memory usage
heap_pattern = re.compile(r'mem_heap_B=(\d+)')
extra_pattern = re.compile(r'mem_heap_extra_B=(\d+)')

try:
    with open(filename, 'r') as file:
        data = file.read()
        
        # Find all matches in the data
        heap_matches = heap_pattern.findall(data)
        extra_matches = extra_pattern.findall(data)

        # Calculate the maximum memory usage in megabytes
        for heap, extra in zip(heap_matches, extra_matches):
            heap_B = int(heap)
            extra_B = int(extra)
            total_memory_B = heap_B + extra_B
            total_memory_MB = total_memory_B / (1024 * 1024)  # Convert bytes to megabytes
            if total_memory_MB > max_memory_usage_MB:
                max_memory_usage_MB = total_memory_MB

    # Print the maximum memory usage
    print(f"Maximum memory usage: {max_memory_usage_MB:.2f} MB")

except FileNotFoundError:
    print(f"Error: The file '{filename}' does not exist.")
except Exception as e:
    print(f"An error occurred: {e}")
