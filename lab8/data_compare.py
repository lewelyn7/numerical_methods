import numpy as np

actual_values = []
with open("result1f.txt", "r") as actual:
    for item in actual:
        actual_values.append(np.complex(float(item.split()[0]), float(item.split()[1])))

# print(actual_values)

expected_values = []
with open("expected1.txt", "r") as expected:
    for item in expected:
        expected_values.append(np.complex(float(item.split()[0]), float(item.split()[1])))

expected_values = np.array(expected_values)
actual_values = np.array(actual_values)

print(np.allclose(actual_values, expected_values))