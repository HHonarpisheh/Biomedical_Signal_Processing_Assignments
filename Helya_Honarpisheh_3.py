import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

def mytransform(img, scale, xshift, yshift, alpha):
    
    Ny, Nx = img.shape

    y = np.arange(Ny)
    x = np.arange(Nx)

    interpolator = RegularGridInterpolator((y, x), img, bounds_error=False, fill_value=0)

    angle = np.radians(alpha)  # Convert degrees to radians
    rotation_matrix = np.array([[np.cos(angle), -np.sin(angle)],
                                [np.sin(angle), np.cos(angle)]])

    new_y = np.linspace(0, Ny - 1, Ny)
    new_x = np.linspace(0, Nx - 1, Nx)
    new_X, new_Y = np.meshgrid(new_x, new_y)

    flat_indices = np.vstack((new_Y.flatten(), new_X.flatten()))

    center_y = Ny / 2.0
    center_x = Nx / 2.0
    flat_indices[0, :] -= center_y
    flat_indices[1, :] -= center_x

    flat_indices *= scale

    rotated_indices = rotation_matrix @ flat_indices

    rotated_indices[0, :] += center_y + yshift
    rotated_indices[1, :] += center_x + xshift

    # Interpolate the rotated image using the transformed coordinates
    rotated_img = interpolator(rotated_indices.T).reshape(Ny, Nx)

    return rotated_img

# Example script to demonstrate the function
if __name__ == "__main__":
    img = plt.imread('/Users/helya.honarpisheh/Library/CloudStorage/OneDrive-Personal/Documents/CCNY/Fall 2024/MedicalImaging/flower.jpg')
    img = img[:, :, 0]  # Convert to grayscale by selecting the first channel

    scale = 1.5
    xshift = 5
    yshift = -3
    alpha = 45  # Rotate by 30 degrees
    imgnew = mytransform(img, scale, xshift, yshift, alpha)

    plt.figure(figsize=(8, 4))

    plt.subplot(1, 2, 1)
    plt.imshow(img, cmap='gray')
    plt.title('Original Image')

    plt.subplot(1, 2, 2)
    plt.imshow(imgnew, cmap='gray')
    plt.title(f'Transformed Image\n(scale={scale}, xshift={xshift}, yshift={yshift}, alpha={alpha} degrees)')

    plt.show()
