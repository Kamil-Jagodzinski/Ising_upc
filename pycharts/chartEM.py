import matplotlib.pyplot as plt

# Function to read data from a text file
def read_data(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        data = [float(line.strip()) for line in lines]
    return data

# Read magnetism data from file
magnetism_data = read_data('avgMagnetism.txt')

# Read energy data from file
energy_data = read_data('energy.txt')

# Plot magnetism data
plt.subplot(2, 1, 1)
plt.plot(magnetism_data, marker='o')
plt.xlabel('Data Point')
plt.ylabel('Magnetism')
plt.title('Magnetism Diagram')

# Save magnetism diagram as PNG file
plt.savefig('magnetism_diagram.png')

# Clear the current figure
plt.clf()

# Plot energy data
plt.subplot(2, 1, 2)
plt.plot(energy_data, marker='o')
plt.xlabel('Data Point')
plt.ylabel('Energy')
plt.title('Energy Diagram')

# Save energy diagram as PNG file
plt.savefig('energy_diagram.png')

# Clear the current figure
plt.clf()