from Image import Image


# example of a catalog containing information about 3 stars
catalog = [[25.65, 31.54, 6.0], [5.979, 14.974, 4.7], [10.77, 17.87, 5.2]]

# initialize the class Image with 50 pixels on each side
myimage = Image(50, 50)

# place the stars of the catalog
myimage.placeStar(catalog)

# place 20 random stars with apparent magnitude between 4 and 5
myimage.place_randomStar(20, 4, 5)

# do a plot of the image with matplotlib.pyplot
myimage.plotImage()

# save the final image into a file in the folder 'output'
myimage.SaveToFile()
