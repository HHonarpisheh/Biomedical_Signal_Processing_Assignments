import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
###
## Slide 10 Figure
# Load the image from the uploaded file path
image_path1 = '/Users/helya.honarpisheh/Library/CloudStorage/OneDrive-Personal/Documents/CCNY/Fall 2024/MedicalImaging/HW1/cell.tif'
cell_img = Image.open(image_path1)
cell_img = np.asarray(cell_img)

# List of resolutions to resize the image to (width x height)
resolutions = [(191, 159), (95, 79), (47, 39), (23, 19)]

# Plot the images at different resolutions
fig1 = plt.figure(figsize=(8, 4))

for i, res in enumerate(resolutions, 1):
    # Resize the image
    resized_img = Image.fromarray(cell_img)
    resized_img = resized_img.resize(res, Image.Resampling.BICUBIC)
    resized_img = np.asarray(resized_img)
    
    # Display the image in a subplot
    plt.subplot(2, 2, i)
    plt.imshow(resized_img, cmap='gray', interpolation='none')
    plt.title(f'{res[1]}x{res[0]}')

plt.tight_layout()

###
from scipy.io import loadmat

## Slide 14 Figure
# Load the image
image_path2 = '/Users/helya.honarpisheh/Library/CloudStorage/OneDrive-Personal/Documents/CCNY/Fall 2024/MedicalImaging/HW1/penny.mat'
penny_img = loadmat(image_path2)['P'] 

# List of intensity levels to display
gray_levels = [2, 4, 8, 16, 32, 256]

fig2 = plt.figure(figsize=(8, 4))

for i, levels in enumerate(gray_levels, 1):
    # Normalize the image to range [0, 1]
    img_normalized = penny_img / 255.0 # chatGPT
    
    # Quantize the image to the desired number of gray levels
    img_quantized = np.round(img_normalized * (levels - 1)) / (levels - 1)
    
    # Rescale back to 0-255 for display
    img_rescaled = (img_quantized * 255).astype(np.uint8) # chatGPT
    
    # Display the quantized image
    plt.subplot(2, 3, i)
    plt.imshow(img_rescaled, cmap='copper', interpolation='none')
    plt.title(f'{levels} gray levels')
    plt.axis('off')

plt.tight_layout()

###
## Slide 19 Figure
# Load the original image
image_path3 = '/Users/helya.honarpisheh/Library/CloudStorage/OneDrive-Personal/Documents/CCNY/Fall 2024/MedicalImaging/HW1/blood1.tif'
blood_img = Image.open(image_path3)
blood_img = np.asarray(blood_img)
blood_img = np.clip(blood_img, 0, 255).astype(np.float64) # chatGPT

# Define the SNR values in dB
snr_values = [np.inf, 9.5, 3.4, -0.12]
fig3, axes = plt.subplots(1, 4, figsize=(8, 4))

# Add noise and plot images with different SNR values
for ax, snr_db in zip(axes, snr_values):
    # Calculate noise variance based on SNR
    snr_linear = 10 ** (snr_db / 10)
    signal_power = np.mean(blood_img ** 2) # chatGPT
    noise_power = signal_power / snr_linear
    
    # Generate Gaussian noise
    noise_std = np.sqrt(noise_power)
    noise = noise_std * np.random.randn(*blood_img.shape)
    
    # Add noise to the image
    noisy_img = blood_img + noise
    
    # Clip pixel values to be in the valid range
    noisy_img = np.clip(noisy_img, 0, 255) # chatGPT
    
    # Plot the noisy image
    ax.imshow(noisy_img, cmap='gray', interpolation='none')
    ax.set_title(f'{snr_db:.2f} dB')
    ax.axis('off')

plt.tight_layout()
###
plt.show()
##
#