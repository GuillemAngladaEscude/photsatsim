from StampLibrary import StampLibrary


# number of partitions of the central pixel fo the stamps
# (the library will have npar_x * npar_y stamps)
npar_x = 10
npar_y = 10

# number of pixels on the stamp
pix_x = 5
pix_y = 5

# desired function to calculate the stamps (see Function module documentation)
func = 0  # 0 -> normalized gaussian

# initialize the library
myLibrary = StampLibrary(func, npar_x, npar_y, pix_x, pix_y)

# generate the library
myLibrary.generateLibrary()

# save the library into a file in the folder 'stamp_libraries'
myLibrary.SaveToFile()

# do a plot of all the stamps of the library
myLibrary.DrawLibrary()

# myLibrary.getError()
