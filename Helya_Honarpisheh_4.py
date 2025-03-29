import numpy as np
import matplotlib.pyplot as plt
from skimage import io, filters, measure, morphology
from scipy.optimize import fmin

# Function to create a plane based on Prabala parameters
def myplane(p, dims):
    [X, Y] = np.meshgrid(np.linspace(0, 1, dims[0]), np.linspace(0, 1, dims[1]))
    return p[0] + p[1] * X + p[2] * Y + p[3] * (X**2) + p[4] * (Y**2)

# Function to fit a plane to the image using Prabala
def myplanefit(Z):
    myerror = lambda p, Z: np.sum((Z - myplane(p, np.shape(Z))) ** 2)
    # Initial parameters set to mean of the image and zeros for linear coefficients
    optimal_params = fmin(myerror, [np.mean(Z), 0, 0, 0, 0], args=(Z,))
    return myplane(optimal_params, np.shape(Z))

# Load the image
image_path = '/Users/helya.honarpisheh/Library/CloudStorage/OneDrive-Personal/Documents/CCNY/Fall 2024/MedicalImaging/HW4/rice.png' 
img = io.imread(image_path, as_gray=True)
img = img / np.max(img)

# Detrend the image to remove background variations
detrended_img = img - myplanefit(img)

# Apply Otsu's thresholding to create a binary image
threshold_value = filters.threshold_otsu(detrended_img)
binary_img = detrended_img > threshold_value

# Create a structuring element for morphological operations
kernel = morphology.square(3)

# Apply morphological operations
eroded_img = morphology.erosion(binary_img, kernel)
dilated_img = morphology.dilation(eroded_img, kernel)
opened_img = morphology.opening(dilated_img, kernel)
closed_img = morphology.closing(opened_img, kernel)

# Label the image and count the number of rice kernels
labelled_img, num_labels = measure.label(closed_img, return_num=True)
properties = measure.regionprops(labelled_img)

# Extract sizes of the detected regions
sizes = [prop.area for prop in properties]

print(f'Total rice kernels detected: {num_labels}')
print(f'Sizes of detected rice kernels: {sizes}')

# Display results
plt.figure(figsize=(12, 6))

plt.subplot(2, 3, 1)
plt.imshow(detrended_img, cmap='gray')
plt.title('Detrended Image')
plt.axis('on')

plt.subplot(2, 3, 2)
plt.imshow(binary_img, cmap='gray')
plt.title('Binary Image')
plt.axis('on')

plt.subplot(2, 3, 3)
plt.imshow(labelled_img, cmap='nipy_spectral')
plt.title('Labelled Image')
plt.axis('on')

plt.subplot(2, 3, 4)
plt.imshow(eroded_img, cmap='gray')
plt.title('Eroded Image')
plt.axis('on')

plt.subplot(2, 3, 5)
plt.imshow(dilated_img, cmap='gray')
plt.title('Dilated Image')
plt.axis('on')

plt.subplot(2, 3, 6)
plt.imshow(opened_img, cmap='gray')
plt.title('Opened Image')
plt.axis('on')

plt.tight_layout()
plt.show()
