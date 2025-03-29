from skimage.io import imread
from skimage.transform import radon, iradon
import numpy as np
import matplotlib.pyplot as plt

image_path = "/Users/helya.honarpisheh/projects/MedImgHWs/CT_image_today.jpg"
initial_image = imread(image_path, as_gray=True)
initial_image = initial_image / np.max(initial_image)  # Normalize to [0, 1]

# Ensure the image is square by padding it to make dimensions equal
size = max(initial_image.shape)
padded_image = np.zeros((size, size))
padded_image[:initial_image.shape[0], :initial_image.shape[1]] = initial_image

# Define projection angles
theta = np.linspace(0., 180., max(padded_image.shape), endpoint=False)

# Actual Radon transform of the padded original image 
g = radon(padded_image, theta=theta, circle=False)

# Initialize image estimate as a uniform image (all ones)
f = np.ones_like(padded_image)

# Number of EM iterations
num_iterations = 50

# EM Iteration loop
for iteration in range(num_iterations):
    # g* estimate of Radon transform
    g_star = radon(f, theta=theta, circle=False)
    
    # Calculate ratio g / g*
    ratio = np.divide(g, g_star + 1e-10, out=np.zeros_like(g), where=g_star != 0)  # Avoid division by zero
    
    # Backproject the ratio with output_size matching f's shape
    backprojection = iradon(ratio, theta=theta, filter_name=None, circle=False, output_size=f.shape[0])
    
    # Update f by multiplying with the backprojection
    f *= backprojection  # Element-wise multiplication to update f


f = f / np.max(f)  # Normalize to [0, 1] range

# Crop f to the original image dimensions
reconstructed_image = f[:initial_image.shape[0], :initial_image.shape[1]]

plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.title("Original Image")
plt.imshow(initial_image, cmap="gray")
plt.axis("on")

plt.subplot(1, 2, 2)
plt.title("Reconstructed Image (after EM)")
plt.imshow(reconstructed_image, cmap="gray")
plt.axis("on")

plt.tight_layout()
plt.show()
