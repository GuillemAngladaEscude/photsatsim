import StampLibrary as sl


# number of stamps = number of partitions
ns_x = 3
ns_y = 3

# number of pixels on the stamp
pix_x = 3
pix_y = 3

myLibrary = sl.StampLibrary(ns_x, ns_x, pix_x, pix_y)
# myLibrary.generateLibrary()
# myLibrary.SaveToFile()
myLibrary.DrawLibrary()
# myLibrary.getError()
