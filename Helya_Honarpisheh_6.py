import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from scipy.signal import convolve2d

def mexican_hat_kernel(size, sigma):
    ax = np.linspace(-(size // 2), size // 2, size)
    xx, yy = np.meshgrid(ax, ax)
    r = np.sqrt(xx**2 + yy**2)  
    kernel = (1 - (r**2 / sigma**2)) * np.exp(-0.5 * (r / sigma)**2)
    return kernel

image_path = '/Users/helya.honarpisheh/Library/CloudStorage/OneDrive-Personal/Documents/CCNY/Fall 2024/MedicalImaging/HW1/cell.tif'
image = Image.open(image_path).convert('L')
image = np.asarray(image)

kernel_size = 21 
sigma = 2 

# Generate the Mexican Hat kernel
mexican_hat = mexican_hat_kernel(kernel_size, sigma)
mexican_hat = mexican_hat / np.sum(mexican_hat) # Normalize the kernel ### chatGPT

# Convolve the image with the Mexican Hat kernel
filtered_image = convolve2d(image, mexican_hat, mode='same', boundary='wrap')

# Convert images to frequency space
F_blurred = np.fft.fft2(filtered_image)
F_filter = np.fft.fft2(mexican_hat, s=filtered_image.shape)

# Set a small threshold value and an offset ### chatGPT
threshold = 1e-10
C = 1e-5
F_filter_with_offset = np.where(np.abs(F_filter) < threshold, C, F_filter)

# Unblur in frequency space with thresholded filter
F_restored = F_blurred / F_filter_with_offset

# Inverse FFT to get the restored image
restored_image = np.abs(np.fft.ifft2(F_restored))
# Normalize the restored image
restored_image = np.clip(restored_image, 0, 255) 

# Display results
plt.figure(figsize=(12, 4))
plt.subplot(1, 3, 1)
plt.title('Original Image')
plt.imshow(image, cmap='gray')

plt.subplot(1, 3, 2)
plt.title('Blurred Image with Mexican Hat')
plt.imshow(filtered_image, cmap='gray')

plt.subplot(1, 3, 3)
plt.title('Restored Image')
plt.imshow(restored_image, cmap='gray')

plt.tight_layout()
plt.show()
