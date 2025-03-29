import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Qt5Agg')  # Ensure smooth animation updates
from skimage.io import imread
from skimage.transform import radon

# Load the brain image (grayscale)
brain_img = imread('/Users/helya.honarpisheh/Library/CloudStorage/OneDrive-Personal/Documents/CCNY/Fall 2024/MedicalImaging/HW5/CT_today_image.jpg', as_gray=True)
rows, cols = brain_img.shape

# Define the projection angles (0 to 179 degrees)
theta = np.linspace(0., 180., max(rows, cols), endpoint=False)

sinogram_manual = np.zeros((len(theta), cols))

# Compute the image center
x_center, y_center = (rows - 1) / 2, (cols - 1) / 2

plt.figure(figsize=(12, 6))

# Display the original image on the left
plt.subplot(1, 3, 1)
plt.imshow(brain_img, cmap='gray')
plt.title('Original Image')
plt.axis('on')

# Activate the subplot for the animated sinogram
plt.subplot(1, 3, 2)
# Loop over every pixel in the image
for i in range(rows):  # For each row
    for j in range(cols):  # For each column
        # Get the pixel value
        pixel_value = brain_img[i, j]

        # Shift coordinates to the image center
        x = i - x_center
        y = j - y_center

        # Loop over all projection angles
        for k, phi in enumerate(theta):
            phi_rad = np.deg2rad(phi)  # Convert to radians

            # Compute the projection position s
            s = x * np.cos(phi_rad) + y * np.sin(phi_rad)

            # Map s to the closest detector bin # chatGPT
            s_idx = int(np.round((s + cols / 2) - 0.5))

            # Ensure the index is within valid bounds # chatGPT
            if 0 <= s_idx < cols:
                # Accumulate the pixel value in the sinogram
                sinogram_manual[k, s_idx] += pixel_value

        # Clear and update the sinogram display
        plt.clf()  # Clear previous plot
        plt.subplot(1, 3, 2)  # Reactivate the subplot
        plt.imshow(sinogram_manual, cmap='gray', aspect='auto')
        plt.title(f'Sinogram Update: Pixel [{i},{j}]')
        plt.xlabel('Projection position (s)')
        plt.ylabel('Projection angle (phi)')

        # Pause to allow matplotlib to render the update
        plt.pause(0.001)
        plt.draw()

# Compute the Radon transform using the built-in function for comparison
sinogram_builtin = radon(brain_img, theta=theta, circle=True)

# Display the built-in Radon output in the third subplot
plt.subplot(1, 3, 3)  
plt.imshow(sinogram_manual, cmap='gray', aspect='auto')
plt.title('Radon() Function Output')
plt.xlabel('Projection position (s)')
plt.ylabel('Projection angle (phi)')

# Final adjustments to the layout and show the figure
plt.tight_layout()
plt.show()
