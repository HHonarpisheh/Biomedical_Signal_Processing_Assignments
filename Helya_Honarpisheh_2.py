import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from scipy.signal import convolve2d


def gaussian_noise(size,sigma):
    ax = np.linspace(-(size // 2), size // 2, size)
    kernel = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * (ax / sigma)**2)
    return kernel

def mexican_hat_kernel(size, sigma):
    # Generate a grid of (x, y) coordinates
    ax = np.linspace(-(size // 2), size // 2, size)
    xx, yy = np.meshgrid(ax, ax)
    r = np.sqrt(xx**2 + yy**2)  # Radial distance from the center *chatgpt
    # Apply the Mexican Hat formula in 2D
    kernel = (1 - (r**2 / sigma**2)) * np.exp(-0.5 * (r / sigma)**2)
    return kernel

def positive_negative_gaussian_kernel(size, sigma):
    # Create a grid of points for the kernel
    ax = np.linspace(-(size // 2), size // 2, size)
    xx, yy = np.meshgrid(ax, ax)
    # Define the 2D positive Gaussian function
    y_positive = (1 / (2 * np.pi * sigma ** 2)) * np.exp(-((xx + size // 4) ** 2 + (yy - 0) ** 2) / (2 * sigma ** 2))
    # Define the 2D negative Gaussian function
    y_negative = -(1 / (2 * np.pi * sigma ** 2)) * np.exp(-((xx - size // 4) ** 2 + (yy - 0) ** 2) / (2 * sigma ** 2))
    # Combine the two to form the kernel
    kernel = y_positive + y_negative
    return kernel


image_path = '/Users/helya.honarpisheh/Library/CloudStorage/OneDrive-Personal/Documents/CCNY/Fall 2024/MedicalImaging/HW1/cell.tif'
image = Image.open(image_path)
image = np.asarray(image)

kernel_size = 40  # Size of the kernel
sigma = 2        # Standard deviation for the Gaussian


# Generate the 1D Gaussian noise
gaussian_kernel = gaussian_noise(kernel_size, sigma)
gaussian_kernel_2d = np.outer(gaussian_kernel, gaussian_kernel)  # Create 2D kernel from 1D *chatgpt

# Convolve the noisy image with the Gaussian kernel
filtered_noisy_image = convolve2d(image, gaussian_kernel_2d, mode='same')


plt.figure(figsize=(8, 4))

plt.subplot(1, 3, 1)
plt.title('Original Image')
plt.imshow(image, cmap='gray')
plt.axis('on')

plt.subplot(1, 3, 2)
plt.title('Gaussian Kernel')
plt.imshow(gaussian_kernel_2d, cmap='gray')
plt.colorbar()
plt.axis('on')

plt.subplot(1, 3, 3)
plt.title('Noisy Image')
plt.imshow(filtered_noisy_image, cmap='gray')
plt.axis('on')

plt.tight_layout()

# Generate the Mexican Hat kernel
mexican_hat = mexican_hat_kernel(kernel_size, sigma)

# Convolve the image with the Mexican Hat kernel
filtered_image = convolve2d(image, mexican_hat, mode='same')


plt.figure(figsize=(8, 4))

plt.subplot(1, 3, 1)
plt.title('Original Image')
plt.imshow(image, cmap='gray')
plt.axis('on')

plt.subplot(1, 3, 2)
plt.title('Mexican Hat Kernel')
plt.imshow(mexican_hat, cmap='gray')
plt.colorbar()
plt.axis('on')

plt.subplot(1, 3, 3)
plt.title('Noisy Image')
plt.imshow(filtered_image, cmap='gray')
plt.axis('on')

plt.tight_layout()

# Generate the positive negative gaussian kernel
kernel = positive_negative_gaussian_kernel(kernel_size, sigma)

# Convolve the image with the kernel
noisy_image = convolve2d(image, kernel, mode='same')


plt.figure(figsize=(8, 4))

plt.subplot(1, 3, 1)
plt.imshow(image, cmap='gray')
plt.title('Original Image')
plt.axis('on')

plt.subplot(1, 3, 2)
plt.imshow(kernel, cmap='gray')
plt.title('+ - Kernel')
plt.colorbar()
plt.axis('on')

plt.subplot(1, 3, 3)
plt.imshow(noisy_image, cmap='gray')
plt.title('Noisy Image')
plt.axis('on')

plt.tight_layout()

plt.show()
