IMAGE_PATH = "/Users/antoinegerardin/Documents/data/rt-cetsa/Data for Nick/20210318 LDHA compound plates/20210318 LDHA compound plate 1 6K cells/1.tif"

from bfio import BioReader, BioWriter
from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np
from skimage.transform import hough_line, hough_line_peaks
from skimage.feature import canny
from skimage import color, img_as_ubyte
from skimage.transform import hough_circle, hough_circle_peaks
from skimage.draw import circle_perimeter
from skimage.transform import hough_ellipse, rotate
from skimage.draw import ellipse_perimeter, disk
from skimage.measure import find_contours, approximate_polygon, subdivide_polygon

from sklearn import preprocessing as pre
import math



def hough_transform_lines(image):
    # Classic straight-line Hough transform
    # Set a precision of 0.5 degree.
    tested_angles = np.linspace(-np.pi / 2, np.pi / 2, 360, endpoint=False)
    h, theta, d = hough_line(image, theta=tested_angles)

    # Generating figure 1
    fig, axes = plt.subplots(1, 3, figsize=(15, 6))
    ax = axes.ravel()

    ax[0].imshow(image, cmap=cm.gray)
    ax[0].set_title('Input image')
    ax[0].set_axis_off()

    angle_step = 0.5 * np.diff(theta).mean()
    d_step = 0.5 * np.diff(d).mean()
    bounds = [
        np.rad2deg(theta[0] - angle_step),
        np.rad2deg(theta[-1] + angle_step),
        d[-1] + d_step,
        d[0] - d_step,
    ]
    ax[1].imshow(np.log(1 + h), extent=bounds, cmap=cm.gray, aspect=1 / 1.5)
    ax[1].set_title('Hough transform')
    ax[1].set_xlabel('Angles (degrees)')
    ax[1].set_ylabel('Distance (pixels)')
    ax[1].axis('image')

    ax[2].imshow(image, cmap=cm.gray)
    ax[2].set_ylim((image.shape[0], 0))
    ax[2].set_axis_off()
    ax[2].set_title('Detected lines')

    for _, angle, dist in zip(*hough_line_peaks(h, theta, d)):
        (x0, y0) = dist * np.array([np.cos(angle), np.sin(angle)])
        ax[2].axline((x0, y0), slope=np.tan(angle + np.pi / 2))

    plt.tight_layout()
    plt.show()

def hough_transform_ellipsis(image):
    edges = canny(image, sigma=2.0, low_threshold=0.55, high_threshold=0.8)

    # Perform a Hough Transform
    # The accuracy corresponds to the bin size of the histogram for minor axis lengths.
    # A higher `accuracy` value will lead to more ellipses being found, at the
    # cost of a lower precision on the minor axis length estimation.
    # A higher `threshold` will lead to less ellipses being found, filtering out those
    # with fewer edge points (as found above by the Canny detector) on their perimeter.
    result = hough_ellipse(edges, accuracy=20, threshold=35000, min_size=8, max_size=20)
    result.sort(order='accumulator')

    # Estimated parameters for the ellipse
    best = list(result[-1])
    yc, xc, a, b = (int(round(x)) for x in best[1:5])
    orientation = best[5]

    # Draw the ellipse on the original image
    cy, cx = ellipse_perimeter(yc, xc, a, b, orientation)
    # Draw the edge (white) and the resulting ellipse (red)
    edges = color.gray2rgb(img_as_ubyte(edges))
    edges[cy, cx] = (250, 0, 0)

    plt.imshow(edges)

    plt.show()

def hough_transform_ovals(image):
    edges = canny(image, sigma=3, low_threshold=10, high_threshold=50)

    # plt.imshow(edges, cmap='gray')
    # plt.show()


    # Detect two radii
    hough_radii = np.arange(4, 20, 1)
    hough_res = hough_circle(edges, hough_radii)

    # Select the most prominent circles
    accums, cx, cy, radii = hough_circle_peaks(hough_res, hough_radii, total_num_peaks=900)

    wells = np.sort(np.column_stack((cx, cy)),axis=0)
    
    # uniques = np.unique(wells)


    min_x = np.argmin(cx)
    max_x = np.argmax(cx)
    min_y = np.argmin(cy)
    max_y = np.argmax(cy)

    corners = [[cx[min_x],cy[min_x]], [cx[max_x],cy[max_x]], [cx[min_y],cy[min_y]], [cx[max_y],cy[max_y]]]

    # corners[0] - corners[1]

    

    # rect = approximate_polygon(corners_arr, tolerance=0.02)
    # plt.imshow(rect, cmap=plt.cm.gray)
    # plt.show()

    corners_top = sorted(corners , key=lambda c : c[1])
    top_left, top_right = corners_top[0], corners_top[1]
    slope = (top_left[1] - top_right[1]) / (top_left[0] - top_right[0]);
    angle = math.degrees(math.atan(slope));

    corners_arr = np.array(corners)
    center_x, center_y = np.sum(corners_arr, axis=0) / 4

    print(f"rotating image by {angle} degrees!")

    # image = rotate(image, 45, center=(center_x, center_y))
   # image = rotate(image, angle, center=top_left)
    image = rotate(image, angle)

    plt.imshow(image, cmap=plt.cm.gray)
    plt.show()

    

    print(f"""
        corners: {min_x}({cx[min_x],cy[min_x]}),
        {max_x}({cx[max_x],cy[max_x]}),
        {min_y}({cx[max_x],cy[max_x]}), {max_y}({cx[max_y],cy[max_y]})
        """)

    # Draw them
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10, 4))
    
    image_norm = (image-np.min(image))/(np.max(image)-np.min(image))
    # image_norm = pre.MinMaxScaler().fit_transform(image)

    id = 1

    plt.imshow(image_norm, cmap=plt.cm.gray)
    plt.show()

    mask = np.zeros(image.shape)
    for center_y, center_x, radius in zip(cy, cx, radii):
        rr, cc = disk((center_y, center_x), radius, shape=image.shape)
        mask[rr, cc] = id
        # id += 1

    # mask = rotate(mask, angle, center=top_left)
    mask = rotate(mask, angle)
    ax.imshow(mask, cmap=plt.cm.gray)
    plt.show()

    # filtered_image = np.where(image, mask, 0)
    # filtered_image = image_norm[mask.astype(bool)].reshape(image.shape)
    filtered_image = np.ma.masked_where(mask == 1, image_norm)
    plt.imshow(filtered_image, cmap=plt.cm.gray)
    plt.show()

    unqID,idx,IDsums = np.unique(mask,return_counts=True,return_inverse=True)
    value_sums = np.bincount(idx,image_norm.ravel())
    {i:(IDsums[itr],value_sums[itr]) for itr,i in enumerate(unqID)}

    image = color.gray2rgb(image_norm)
    for center_y, center_x, radius in zip(cy, cx, radii):
        circy, circx = circle_perimeter(center_y, center_x, radius, shape=image.shape)
        image[circy, circx] = (220, 20, 20)

    ax.imshow(image, cmap=plt.cm.gray)
    plt.show()


def process_plate_method1():
    with BioReader(IMAGE_PATH) as br:
        print(br.shape)
        print(br._DIMS)
        bpp = br.bpp
        data = br[:].squeeze()
        print(data[0:200,0:200])
        print(data.shape)
        print(data.ndim)
        print(f" max : {2 **( bpp * 8)}")
        # plt.imshow(data, interpolation='nearest')
        #plt.imshow(data /  2 **( bpp * 8), cmap='gray')
        # equivalent imshow can rescale no problem
        plt.imshow(data, cmap='gray')
        # plt.show()
        data /  2 **( bpp * 8)
        # hough_transform_lines(data)
        hough_transform_ovals(data)
        # hough_transform_ellipsis(data)



if __name__ == "__main__":
    process_plate_method1()