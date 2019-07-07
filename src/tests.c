#include "connected_layer.h"
#include "convolutional_layer.h"
#include "maxpool_layer.h"
#include "network.h"
#include "image.h"
#include "parser.h"
#include "data.h"
#include "matrix.h"
#include "utils.h"

#include <time.h>
#include <stdlib.h>
#include <stdio.h>


double error_network(network net, matrix m, double **truth)
{
	int i;
	int correct = 0;
	int k = get_network_output_size(net);
	for (i = 0; i < m.rows; ++i)
	{
		forward_network(net, m.vals[i]);
		double *out = get_network_output(net);
		int guess = max_index(out, k);
		if (truth[i][guess])
			++correct;
	}
	return (double)correct / m.rows;
}


double **one_hot(double *a, int n, int k)
{
    int i;
    double **t = calloc(n, sizeof(double*));
    for(i = 0; i < n; ++i)
	{
        t[i] = calloc(k, sizeof(double));

		// get the true number from array a
        int index = (int)a[i];
        t[i][index] = 1;
    }
    return t;
}

void test_nist()
{
    srand(999999);
    network net = parse_network_cfg("nist.cfg");

	// m.rows means how many numbers
    matrix m = csv_to_matrix("mnist/mnist_train.csv");
    matrix test = csv_to_matrix("mnist/mnist_test.csv");

	// pop the first column, because it's the target number.
    double *truth_1d = pop_column(&m, 0);

	// create a new array to save the first column
    double **truth = one_hot(truth_1d, m.rows, 10);

	// pop the first column
    double *test_truth_1d = pop_column(&test, 0);
    double **test_truth = one_hot(test_truth_1d, test.rows, 10);

    int i,j;
    clock_t start = clock(), end;
    for(i = 0; i < test.rows; ++i){
        normalize_array(test.vals[i], 28*28);
    }
    for(i = 0; i < m.rows; ++i){
        normalize_array(m.vals[i], 28*28);
    }


	// training， 1. 使用训练数据
    int count = 0;
    double lr = .0005;
    while(++count <= 300)
	{
        //lr *= .99;
        int index = 0;
        int correct = 0;
		int number = 1000;
        for(i = 0; i < number; ++i)
		{
			index = rand() % m.rows;

			// 2. 得到预测值,  output = active * (input * weight + bias)
            forward_network(net, m.vals[index]);
            double *out = get_network_output(net);
            double *delta = get_network_delta(net);

            int max_i = 0;
            double max = out[0];
            for(j = 0; j < 10; ++j)
			{
				// get delta， 3， 得到差值
                delta[j] = truth[index][j]-out[j];

				// find the max
                if(out[j] > max)
				{
                    max = out[j];
                    max_i = j;
                }
            }

			// compare result with the true number
            if(truth[index][max_i])
				++correct;

			// learn, loss function， 3， 计算loss， 保存结果到bias_updates
            learn_network(net, m.vals[index]);

			// 4, 获得最终biases和kernel(weight)
            update_network(net, lr);
        }

        //print_network(net);
		//image input = double_to_image(28, 28, 1, m.vals[index]);

        //show_image(input, "Input");
        //image o = get_network_image(net);
		//show_image(o, "Output_Original");
        //show_image_collapsed(o, "Output");

		// show the conv layer?
        //visualize_network(net);

        //cvWaitKey(10);
        //double test_acc = error_network(net, m, truth);
        fprintf(stderr, "\n%5d: %f %f\n\n",count, (double)correct/number, lr);
        if(count % (m.rows/number) == 0)
			lr /= 2; 
    }

	// test
	double train_acc = error_network(net, m, truth);
	fprintf(stderr, "\nTRAIN: %f\n", train_acc);
	double test_acc = error_network(net, test, test_truth);
	fprintf(stderr, "TEST: %f\n\n", test_acc);
	printf("%d, %f, %f\n", count, train_acc, test_acc);
    end = clock();
    //printf("Neural Net Learning: %lf seconds\n", (double)(end-start)/CLOCKS_PER_SEC);
}


int main()
{
    test_nist();
    return 0;
}
