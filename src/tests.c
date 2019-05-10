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

void test_convolve()
{
    image dog = load_image("dog.jpg");
    printf("dog channels %d\n", dog.c);
    image kernel = make_random_image(3,3,dog.c);
    image edge = make_image(dog.h, dog.w, 1);
    int i;
    clock_t start = clock(), end;
    for(i = 0; i < 1000; ++i){
        convolve(dog, kernel, 1, 0, edge, 1);
    }
    end = clock();
    printf("Convolutions: %lf seconds\n", (double)(end-start)/CLOCKS_PER_SEC);
    show_image_layers(edge, "Test Convolve");
}

void test_color()
{
    image dog = load_image("test_color.png");
    show_image_layers(dog, "Test Color");
}

void test_convolutional_layer()
{
    srand(0);
    image dog = load_image("dog.jpg");
    int i;
    int n = 3;
    int stride = 1;
    int size = 3;
    convolutional_layer layer = *make_convolutional_layer(dog.h, dog.w, dog.c, n, size, stride, RELU);
    char buff[256];
    for(i = 0; i < n; ++i) {
        sprintf(buff, "Kernel %d", i);
        show_image(layer.kernels[i], buff);
    }
    forward_convolutional_layer(layer, dog.data);
    
    image output = get_convolutional_image(layer);
    maxpool_layer mlayer = *make_maxpool_layer(output.h, output.w, output.c, 2);
    forward_maxpool_layer(mlayer, layer.output);

    show_image_layers(get_maxpool_image(mlayer), "Test Maxpool Layer");
}

void verify_convolutional_layer()
{
    srand(0);
    int i;
    int n = 1;
    int stride = 1;
    int size = 3;
    double eps = .00000001;
    image test = make_random_image(5,5, 1);
    convolutional_layer layer = *make_convolutional_layer(test.h,test.w,test.c, n, size, stride, RELU);
    image out = get_convolutional_image(layer);
    double **jacobian = calloc(test.h*test.w*test.c, sizeof(double));
    
    forward_convolutional_layer(layer, test.data);
    image base = copy_image(out);

    for(i = 0; i < test.h*test.w*test.c; ++i){
        test.data[i] += eps;
        forward_convolutional_layer(layer, test.data);
        image partial = copy_image(out);
        subtract_image(partial, base);
        scale_image(partial, 1/eps);
        jacobian[i] = partial.data;
        test.data[i] -= eps;
    }
    double **jacobian2 = calloc(out.h*out.w*out.c, sizeof(double));
    image in_delta = make_image(test.h, test.w, test.c);
    image out_delta = get_convolutional_delta(layer);
    for(i = 0; i < out.h*out.w*out.c; ++i){
        out_delta.data[i] = 1;
        backward_convolutional_layer(layer, test.data, in_delta.data);
        image partial = copy_image(in_delta);
        jacobian2[i] = partial.data;
        out_delta.data[i] = 0;
    }
    int j;
    double *j1 = calloc(test.h*test.w*test.c*out.h*out.w*out.c, sizeof(double));
    double *j2 = calloc(test.h*test.w*test.c*out.h*out.w*out.c, sizeof(double));
    for(i = 0; i < test.h*test.w*test.c; ++i){
        for(j =0 ; j < out.h*out.w*out.c; ++j){
            j1[i*out.h*out.w*out.c + j] = jacobian[i][j];
            j2[i*out.h*out.w*out.c + j] = jacobian2[j][i];
            printf("%f %f\n", jacobian[i][j], jacobian2[j][i]);
        }
    }


    image mj1 = double_to_image(test.w*test.h*test.c, out.w*out.h*out.c, 1, j1);
    image mj2 = double_to_image(test.w*test.h*test.c, out.w*out.h*out.c, 1, j2);
    printf("%f %f\n", avg_image_layer(mj1,0), avg_image_layer(mj2,0));
    show_image(mj1, "forward jacobian");
    show_image(mj2, "backward jacobian");
    
}

void test_load()
{
    image dog = load_image("dog.jpg");
    show_image(dog, "Test Load");
    show_image_layers(dog, "Test Load");
}
void test_upsample()
{
    image dog = load_image("dog.jpg");
    int n = 3;
    image up = make_image(n*dog.h, n*dog.w, dog.c);
    upsample_image(dog, n, up);
    show_image(up, "Test Upsample");
    show_image_layers(up, "Test Upsample");
}

void test_rotate()
{
    int i;
    image dog = load_image("dog.jpg");
    clock_t start = clock(), end;
    for(i = 0; i < 1001; ++i){
        rotate_image(dog);
    }
    end = clock();
    printf("Rotations: %lf seconds\n", (double)(end-start)/CLOCKS_PER_SEC);
    show_image(dog, "Test Rotate");

    image random = make_random_image(3,3,3);
    show_image(random, "Test Rotate Random");
    rotate_image(random);
    show_image(random, "Test Rotate Random");
    rotate_image(random);
    show_image(random, "Test Rotate Random");
}

void test_parser()
{
    network net = parse_network_cfg("test_parser.cfg");
    double input[1];
    int count = 0;
        
    double avgerr = 0;
    while(++count < 100000000){
        double v = ((double)rand()/RAND_MAX);
        double truth = v*v;
        input[0] = v;
        forward_network(net, input);
        double *out = get_network_output(net);
        double *delta = get_network_delta(net);
        double err = pow((out[0]-truth),2.);
        avgerr = .99 * avgerr + .01 * err;
        if(count % 1000000 == 0) printf("%f %f :%f AVG %f \n", truth, out[0], err, avgerr);
        delta[0] = truth - out[0];
        learn_network(net, input);
        update_network(net, .001);
    }
}

void test_data()
{
    char *labels[] = {"cat","dog"};
    batch train = random_batch("train_paths.txt", 101,labels, 2);
    show_image(train.images[0], "Test Data Loading");
    show_image(train.images[100], "Test Data Loading");
    show_image(train.images[10], "Test Data Loading");
    free_batch(train);
}

void test_full()
{
    network net = parse_network_cfg("full.cfg");
    srand(0);
    int i = 0;
    char *labels[] = {"cat","dog"};
    while(i++ < 1000 || 1){
        batch train = random_batch("train_paths.txt", 1000, labels, 2);
        train_network_batch(net, train);
        free_batch(train);
        printf("Round %d\n", i);
    }
}

double error_network(network net, matrix m, double **truth)
{
    int i;
    int correct = 0;
    int k = get_network_output_size(net);
    for(i = 0; i < m.rows; ++i){
        forward_network(net, m.vals[i]);
        double *out = get_network_output(net);
        int guess = max_index(out, k);
        if(truth[i][guess]) ++correct;
    }
    return (double)correct/m.rows;
}

double **one_hot(double *a, int n, int k)
{
    int i;
    double **t = calloc(n, sizeof(double*));
    for(i = 0; i < n; ++i){
        t[i] = calloc(k, sizeof(double));
        int index = (int)a[i];
        t[i][index] = 1;
    }
    return t;
}

void test_nist()
{
    srand(999999);
    network net = parse_network_cfg("nist.cfg");
    matrix m = csv_to_matrix("mnist/mnist_train_100.csv");
    matrix test = csv_to_matrix("mnist/mnist_test_10.csv");
    double *truth_1d = pop_column(&m, 0);
    double **truth = one_hot(truth_1d, m.rows, 10);
    double *test_truth_1d = pop_column(&test, 0);
    double **test_truth = one_hot(test_truth_1d, test.rows, 10);
    int i,j;
    clock_t start = clock(), end;
    for(i = 0; i < test.rows; ++i){
        normalize_array(test.vals[i], 28*28);
        //scale_array(m.vals[i], 28*28, 1./255.);
        //translate_array(m.vals[i], 28*28, -.1);
    }
    for(i = 0; i < m.rows; ++i){
        normalize_array(m.vals[i], 28*28);
        //scale_array(m.vals[i], 28*28, 1./255.);
        //translate_array(m.vals[i], 28*28, -.1);
    }
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
            index = rand()%m.rows;
            forward_network(net, m.vals[index]);
            double *out = get_network_output(net);
            double *delta = get_network_delta(net);
            int max_i = 0;
            double max = out[0];
            for(j = 0; j < 10; ++j)
			{
                delta[j] = truth[index][j]-out[j];
                if(out[j] > max){
                    max = out[j];
                    max_i = j;
                }
            }
            if(truth[index][max_i]) ++correct;
            learn_network(net, m.vals[index]);
            update_network(net, lr);
        }

        print_network(net);
        image input = double_to_image(28,28,1, m.vals[index]);
        //show_image(input, "Input");
        image o = get_network_image(net);
        //show_image_collapsed(o, "Output");
        visualize_network(net);
        //cvWaitKey(10);
        //double test_acc = error_network(net, m, truth);
        fprintf(stderr, "\n%5d: %f %f\n\n",count, (double)correct/number, lr);
        if(count % 10 == 0 && 0)
		{
            double train_acc = error_network(net, m, truth);
            fprintf(stderr, "\nTRAIN: %f\n", train_acc);
            double test_acc = error_network(net, test, test_truth);
            fprintf(stderr, "TEST: %f\n\n", test_acc);
            printf("%d, %f, %f\n", count, train_acc, test_acc);
        }

        if(count % (m.rows/number) == 0) lr /= 2; 
    }


            double train_acc = error_network(net, m, truth);
            fprintf(stderr, "\nTRAIN: %f\n", train_acc);
            double test_acc = error_network(net, test, test_truth);
            fprintf(stderr, "TEST: %f\n\n", test_acc);
            printf("%d, %f, %f\n", count, train_acc, test_acc);
    end = clock();
    //printf("Neural Net Learning: %lf seconds\n", (double)(end-start)/CLOCKS_PER_SEC);
}

void test_kernel_update()
{
    srand(0);
    double delta[] = {.1};
    double input[] = {.3, .5, .3, .5, .5, .5, .5, .0, .5};
    double kernel[] = {1,2,3,4,5,6,7,8,9};
    convolutional_layer layer = *make_convolutional_layer(3, 3, 1, 1, 3, 1, LINEAR);
    layer.kernels[0].data = kernel;
    layer.delta = delta;
    learn_convolutional_layer(layer, input);
    print_image(layer.kernels[0]);
    print_image(get_convolutional_delta(layer));
    print_image(layer.kernel_updates[0]);

}

void test_random_classify()
{
    network net = parse_network_cfg("connected.cfg");
    matrix m = csv_to_matrix("train.csv");
    matrix ho = hold_out_matrix(&m, 2500);
    double *truth = pop_column(&m, 0);
    double *ho_truth = pop_column(&ho, 0);
    int i;
    clock_t start = clock(), end;
    int count = 0;
    while(++count <= 300){
        for(i = 0; i < m.rows; ++i){
            int index = rand()%m.rows;
            //image p = double_to_image(1690,1,1,m.vals[index]);
            //normalize_image(p);
            forward_network(net, m.vals[index]);
            double *out = get_network_output(net);
            double *delta = get_network_delta(net);
            //printf("%f\n", out[0]);
            delta[0] = truth[index] - out[0];
            // printf("%f\n", delta[0]);
            //printf("%f %f\n", truth[index], out[0]);
            learn_network(net, m.vals[index]);
            update_network(net, .00001);
        }
        //double test_acc = error_network(net, m, truth);
        //double valid_acc = error_network(net, ho, ho_truth);
        //printf("%f, %f\n", test_acc, valid_acc);
        //fprintf(stderr, "%5d: %f Valid: %f\n",count, test_acc, valid_acc);
        //if(valid_acc > .70) break;
    }
    end = clock();
    FILE *fp = fopen("submission/out.txt", "w");
    matrix test = csv_to_matrix("test.csv");
    truth = pop_column(&test, 0);
    for(i = 0; i < test.rows; ++i){
        forward_network(net, test.vals[i]);
        double *out = get_network_output(net);
        if(fabs(out[0]) < .5) fprintf(fp, "0\n");
        else fprintf(fp, "1\n");
    }
    fclose(fp);
    printf("Neural Net Learning: %lf seconds\n", (double)(end-start)/CLOCKS_PER_SEC);
}

void test_random_preprocess()
{
    FILE *file = fopen("train.csv", "w");
    char *labels[] = {"cat","dog"};
    int i,j,k;
    srand(0);
    network net = parse_network_cfg("convolutional.cfg");
    for(i = 0; i < 100; ++i){
        printf("%d\n", i);
        batch part = get_batch("train_paths.txt", i, 100, labels, 2);
        for(j = 0; j < part.n; ++j){
            forward_network(net, part.images[j].data);
            double *out = get_network_output(net);
            fprintf(file, "%f", part.truth[j][0]);
            for(k = 0; k < get_network_output_size(net); ++k){
                fprintf(file, ",%f", out[k]);
            }
            fprintf(file, "\n");
        }
        free_batch(part);
    }
}

int main()
{
    //test_kernel_update();
    test_nist();
    //test_full();
    //test_random_preprocess();
    //test_random_classify();
    //test_parser();
    //test_backpropagate();
    //test_ann();
    //test_convolve();
    //test_upsample();
    //test_rotate();
    //test_load();
    //test_network();
    //test_convolutional_layer();
    //verify_convolutional_layer();
    //test_color();
    //cvWaitKey(0);
    return 0;
}
