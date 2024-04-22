import matplotlib.pyplot as plt
from matplotlib import transforms

IMAGE_PATH = "/Users/antoinegerardin/Documents/data/rt-cetsa/Data for Nick/20210318 LDHA compound plates/20210318 LDHA compound plate 1 6K cells/1.tif"


img = plt.imread(IMAGE_PATH)

fig = plt.figure()
ax = fig.add_subplot(111)

rotation_in_degrees = -20

tr = transforms.Affine2D().rotate_deg(rotation_in_degrees)

ax.imshow(img, transform=tr + ax.transData)
plt.show()