import numpy as np
from scipy.optimize import fmin
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from skimage import io

# Define the paths for the images
image_1_path = '/Users/helya.honarpisheh/Library/CloudStorage/OneDrive-Personal/Documents/CCNY/Fall 2024/MedicalImaging/HW3/head-1.png'
image_2_path = '/Users/helya.honarpisheh/Library/CloudStorage/OneDrive-Personal/Documents/CCNY/Fall 2024/MedicalImaging/HW3/head-2.png'

image_1 = io.imread(image_1_path, as_gray=True)
image_2 = io.imread(image_2_path, as_gray=True)

def mytransform(img, scale, xshift, yshift, alpha):
    Ny, Nx = img.shape
    
    y = np.arange(Ny)
    x = np.arange(Nx)
    
    interpolator = RegularGridInterpolator((y, x), img, bounds_error=False, fill_value=0) # chatGPT
    
    # Convert degrees to radians for rotation
    angle = np.radians(alpha)
    rotation_matrix = np.array([[np.cos(angle), -np.sin(angle)],
                                [np.sin(angle), np.cos(angle)]])
    
    new_y = np.linspace(0, Ny - 1, Ny)
    new_x = np.linspace(0, Nx - 1, Nx)
    new_X, new_Y = np.meshgrid(new_x, new_y)
    
    flat_indices = np.vstack((new_Y.flatten(), new_X.flatten())) # chatGPT
    
    # Center the coordinates for rotation
    center_y = Ny / 2.0
    center_x = Nx / 2.0
    flat_indices[0, :] -= center_y
    flat_indices[1, :] -= center_x
    
    # Apply scaling
    flat_indices *= scale
    
    # Apply rotation
    rotated_indices = rotation_matrix @ flat_indices
    
    # Apply shifting
    rotated_indices[0, :] += center_y + yshift
    rotated_indices[1, :] += center_x + xshift
    
    img_new = interpolator((rotated_indices[0, :], rotated_indices[1, :])).reshape(Ny, Nx)
    
    return img_new

def error_function(params, img_fixed, img_moving):
    scale, xshift, yshift, alpha = params
    img_transformed = mytransform(img_moving, scale, xshift, yshift, alpha)
    # Compute mean squared error
    mse = np.mean((img_fixed - img_transformed) ** 2)
    return mse

def coregister_images(img_fixed, img_moving):
    num_restarts = 10
    best_result = None
    best_mse = float('inf')

    for _ in range(num_restarts):
        # Generate random initial parameters
        initial_params = [np.random.uniform(0.5, 1.5),  # Scale
                        np.random.uniform(-10, 10),  # X Shift
                        np.random.uniform(-10, 10),  # Y Shift
                        np.random.uniform(-45, 45)]  # Alpha 

        # Use fmin instead of minimize
        result = fmin(error_function, initial_params, args=(img_fixed, img_moving), disp=False)

        current_mse = error_function(result, img_fixed, img_moving)
        if current_mse < best_mse:
            best_mse = current_mse
            best_result = result
    
    # Use best_result for the optimal parameters
    img_registered = mytransform(img_moving, *best_result)
    
    return img_registered, best_result

# Perform coregistration of the two images
registered_image, optimal_params = coregister_images(image_1, image_2)

optimal_params_dict = {
    'Optimal Scale': optimal_params[0],
    'Optimal X Shift': optimal_params[1],
    'Optimal Y Shift': optimal_params[2],
    'Optimal Rotation Angle': optimal_params[3]
}

print(optimal_params_dict)

plt.figure(figsize=(12, 4))

plt.subplot(1, 3, 1)
plt.imshow(image_1, cmap='gray')
plt.title('Fixed Image')
plt.axis('on')

plt.subplot(1, 3, 2)
plt.imshow(registered_image, cmap='gray')
plt.title('Registered Moving Image')
plt.axis('on')

plt.subplot(1, 3, 3)
plt.imshow(np.abs(image_1 - registered_image), cmap='gray')
plt.title('Difference Image')
plt.axis('on')

plt.show()
