from Image import Image


positions_stars = [[25.65, 31.54, 7.0], [5.979, 14.974, 11.7], [10.77, 17.87, 10.2]]
# positions_stars = [[2.96, 1.97, 10000.0]]
# mylibrary = np.load('Library_7x7.npy')

myimage = Image(50, 50)
# myimage.initLibraryWithStampSize(3)
# myimage.setPixelsStamp(7,7)
# myimage.setLibrary(mylibrary)
# myimage.placeStar(positions_stars)
myimage.place_randomStar(20, 4, 5)
# myimage.test(135)

myimage.plotImage()
# myimage.SaveToFile()



