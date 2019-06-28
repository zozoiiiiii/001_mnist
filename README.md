# what is it
	train numbers which writed by hand, then test

# how it works

#### net structure
![Layers](img/layers.png)


#### one hot
![One Hot](img/one_hot.png)

# reference

#### yolo
 download code(2013-12-06): https://github.com/pjreddie/darknet

#### mnist dataset
1. original data:
	train images: train-images-idx3-ubyte.gz (60000)
	train labels: train-labels-idx1-ubyte.gz (60000)
	test images: t10k-images-idx3-ubyte.gz (10000)
	test labels: t10k-labels-idx1-ubyte.gz (10000)

2. yolo data: https://pjreddie.com/projects/mnist-in-csv/
	mnist_train.csv
	mnist_test.csv
	
	format: everyline in the csv file is one object:  label + image

#### reference
    https://www.cnblogs.com/bestExpert/p/9291185.html
	https://github.com/makeyourownneuralnetwork/makeyourownneuralnetwork
	