import numpy as np
import matplotlib.pyplot as plt
from skimage.io import imread
from skimage.transform import radon
from scipy.interpolate import interp1d

brain_img = imread('/Users/helya.honarpisheh/Library/CloudStorage/OneDrive-Personal/Documents/CCNY/Fall 2024/MedicalImaging/HW5/CT_today_image.jpg', as_gray=True)
I, J = brain_img.shape

phi = np.linspace(0, 180, max(I, J), endpoint=False)
sinogram = radon(brain_img, theta=phi, circle=False)

G = np.fft.fft(sinogram, axis=0)
k = G.shape[0]
H = np.abs(np.fft.fftfreq(k) * k).reshape(-1, 1) * np.ones((1, sinogram.shape[1]))
filtered_sinogram = np.real(np.fft.ifft(G * H, axis=0))

b = np.zeros((I, J))

x, y = np.meshgrid(np.arange(J) - J / 2, np.arange(I) - I / 2)

for i, angle in enumerate(phi):
    # Compute shifted coordinates for each projection angle
    angle_rad = np.deg2rad(angle)
    s = x * np.cos(angle_rad) - y * np.sin(angle_rad)
    s += filtered_sinogram.shape[0] / 2  # Shift to match sinogram index

    # Interpolate for each projection line in the sinogram
    interp = interp1d(np.arange(filtered_sinogram.shape[0]), filtered_sinogram[:, i], bounds_error=False, fill_value=0)
    b += interp(s.ravel()).reshape(I, J)

# Normalize the backprojected image ###suggested by chatGPT
b -= np.min(b)  # Set minimum to 0
b /= np.max(b)  # Scale to [0, 1] for display

plt.figure(figsize=(16, 4))

# Original Image
plt.subplot(1, 4, 1)
plt.imshow(brain_img, cmap='gray')
plt.title('Original Image')
plt.axis('on')

# Sinogram
plt.subplot(1, 4, 2)
plt.imshow(sinogram, cmap='gray', extent=(phi.min(), phi.max(), -J / 2, J / 2))
plt.title('Sinogram')
plt.xlabel('Projection angle (phi)')
plt.ylabel('Projection position (s)')

# Filtered Sinogram
plt.subplot(1, 4, 3)
plt.imshow(filtered_sinogram, cmap='gray', extent=(phi.min(), phi.max(), -J / 2, J / 2))
plt.title('Filtered Sinogram')
plt.xlabel('Projection angle (phi)')
plt.ylabel('Projection position (s)')

# Reconstructed Image using Back-projection
plt.subplot(1, 4, 4)
plt.imshow(b, cmap='gray')
plt.title('Reconstructed Image')
plt.axis('on')

plt.tight_layout()
plt.show()
