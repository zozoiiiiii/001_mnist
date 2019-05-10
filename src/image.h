#ifndef IMAGE_H
#define IMAGE_H

//#include "opencv2/highgui/highgui_c.h"
//#include "opencv2/imgproc/imgproc_c.h"
typedef struct {
    int h;
    int w;
    int c;
    double *data;
} image;

void scale_image(image m, double s);
void add_scalar_image(image m, double s);
void normalize_image(image p);
void z_normalize_image(image p);
void threshold_image(image p, double t);
void zero_image(image m);
void rotate_image(image m);
void subtract_image(image a, image b);
double avg_image_layer(image m, int l);
void embed_image(image source, image dest, int h, int w);
image collapse_image_layers(image source, int border);

void show_image(image p, char *name);
void show_image_layers(image p, char *name);
void show_image_collapsed(image p, char *name);
void print_image(image m);

image make_image(int h, int w, int c);
image make_empty_image(int h, int w, int c);
image make_random_image(int h, int w, int c);
image make_random_kernel(int size, int c, double scale);
image double_to_image(int h, int w, int c, double *data);
image copy_image(image p);
image load_image(char *filename);

double get_pixel(image m, int x, int y, int c);
double get_pixel_extend(image m, int x, int y, int c);
void set_pixel(image m, int x, int y, int c, double val);


image get_image_layer(image m, int l);

void two_d_convolve(image m, int mc, image kernel, int kc, int stride, image out, int oc, int edge);
void upsample_image(image m, int stride, image out);
void convolve(image m, image kernel, int stride, int channel, image out, int edge);
void back_convolve(image m, image kernel, int stride, int channel, image out, int edge);
void kernel_update(image m, image update, int stride, int channel, image out, int edge);

void free_image(image m);
#endif

