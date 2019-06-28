#include "image.h"
#include "utils.h"
#include <stdio.h>


#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

int windows = 0;

void subtract_image(image a, image b)
{
    int i;
    for(i = 0; i < a.h*a.w*a.c; ++i) a.data[i] -= b.data[i];
}

void embed_image(image source, image dest, int h, int w)
{
    int i,j,k;
    for(k = 0; k < source.c; ++k){
        for(i = 0; i < source.h; ++i){
            for(j = 0; j < source.w; ++j){
                double val = get_pixel(source, i,j,k);
                set_pixel(dest, h+i, w+j, k, val);
            }
        }
    }
}

image collapse_image_layers(image source, int border)
{
    int h = source.h;
    h = (h+border)*source.c - border;
    image dest = make_image(h, source.w, 1);
    int i;
    for(i = 0; i < source.c; ++i){
        image layer = get_image_layer(source, i);
        int h_offset = i*(source.h+border);
        embed_image(layer, dest, h_offset, 0);
        free_image(layer);
    }
    return dest;
}

void z_normalize_image(image p)
{
    normalize_array(p.data, p.h*p.w*p.c);
}

void normalize_image(image p)
{
    double *min = calloc(p.c, sizeof(double));
    double *max = calloc(p.c, sizeof(double));
    int i,j;
    for(i = 0; i < p.c; ++i) min[i] = max[i] = p.data[i*p.h*p.w];

    for(j = 0; j < p.c; ++j){
        for(i = 0; i < p.h*p.w; ++i){
            double v = p.data[i+j*p.h*p.w];
            if(v < min[j]) min[j] = v;
            if(v > max[j]) max[j] = v;
        }
    }
    for(i = 0; i < p.c; ++i){
        if(max[i] - min[i] < .000000001){
            min[i] = 0;
            max[i] = 1;
        }
    }
    for(j = 0; j < p.c; ++j){
        for(i = 0; i < p.w*p.h; ++i){
            p.data[i+j*p.h*p.w] = (p.data[i+j*p.h*p.w] - min[j])/(max[j]-min[j]);
        }
    }
    free(min);
    free(max);
}

double avg_image_layer(image m, int l)
{
    int i;
    double sum = 0;
    for(i = 0; i < m.h*m.w; ++i){
        sum += m.data[l*m.h*m.w + i];
    }
    return sum/(m.h*m.w);
}

void threshold_image(image p, double t)
{
    int i;
    for(i = 0; i < p.w*p.h*p.c; ++i){
        if(p.data[i] < t) p.data[i] = 0;
    }
}

image copy_image(image p)
{
    image copy = p;
    copy.data = calloc(p.h*p.w*p.c, sizeof(double));
    memcpy(copy.data, p.data, p.h*p.w*p.c*sizeof(double));
    return copy;
}

void save_image_png(image im, const char *name)
{
	char buff[256];
	sprintf(buff, "%s.png", name);
	unsigned char *data = (unsigned char*)calloc(im.w*im.h*im.c, sizeof(char));
	int i, k;
	for (k = 0; k < im.c; ++k) {
		for (i = 0; i < im.w*im.h; ++i) {
			data[i*im.c + k] = (unsigned char)(255 * im.data[i + k * im.w*im.h]);
		}
	}
	int success = stbi_write_png(buff, im.w, im.h, im.c, data, im.w*im.c);
	free(data);
	if (!success) fprintf(stderr, "Failed to write ImageData %s\n", buff);
}

void show_image(image p, char *name)
{
#ifdef OPENCV
	show_image_cv(p, name);
#else
	fprintf(stderr, "Not compiled with OpenCV, saving to %s.png instead\n", name);
	save_image_png(p, name);
#endif
}

void show_image_layers(image p, char *name)
{
    int i;
    char buff[256];
    for(i = 0; i < p.c; ++i){
        sprintf(buff, "%s - Layer %d", name, i);
        image layer = get_image_layer(p, i);
        show_image(layer, buff);
        free_image(layer);
    }
}

void show_image_collapsed(image p, char *name)
{
    image c = collapse_image_layers(p, 1);
    show_image(c, name);
    free_image(c);
}

image make_empty_image(int h, int w, int c)
{
    image out;
    out.h = h;
    out.w = w;
    out.c = c;
    return out;
}

image make_image(int h, int w, int c)
{
    image out = make_empty_image(h,w,c);
    out.data = calloc(h*w*c, sizeof(double));
    return out;
}

image double_to_image(int h, int w, int c, double *data)
{
    image out = make_empty_image(h,w,c);
    out.data = data;
    return out;
}

void zero_image(image m)
{
    memset(m.data, 0, m.h*m.w*m.c*sizeof(double));
}

void zero_channel(image m, int c)
{
    memset(&(m.data[c*m.h*m.w]), 0, m.h*m.w*sizeof(double));
}

void rotate_image(image m)
{
    int i,j;
    for(j = 0; j < m.c; ++j){
        for(i = 0; i < m.h*m.w/2; ++i){
            double swap = m.data[j*m.h*m.w + i];
            m.data[j*m.h*m.w + i] = m.data[j*m.h*m.w + (m.h*m.w-1 - i)];
            m.data[j*m.h*m.w + (m.h*m.w-1 - i)] = swap;
        }
    }
}

image make_random_image(int h, int w, int c)
{
    image out = make_image(h,w,c);
    int i;
    for(i = 0; i < h*w*c; ++i){
        out.data[i] = rand_normal();
    }
    return out;
}

void add_scalar_image(image m, double s)
{
    int i;
    for(i = 0; i < m.h*m.w*m.c; ++i) m.data[i] += s;
}

void scale_image(image m, double s)
{
    int i;
    for(i = 0; i < m.h*m.w*m.c; ++i) m.data[i] *= s;
}

image make_random_kernel(int size, int c, double scale)
{
    int pad;
    if((pad=(size%2==0))) ++size;
    image out = make_random_image(size,size,c);
    scale_image(out, scale);
    int i,k;
    if(pad){
        for(k = 0; k < out.c; ++k){
            for(i = 0; i < size; ++i) {
                set_pixel(out, i, 0, k, 0);
                set_pixel(out, 0, i, k, 0);
            }
        }
    }
    return out;
}


image load_image_stb(char *filename, int channels)
{
    int w, h, c;
    unsigned char *data = stbi_load(filename, &w, &h, &c, channels);
    if (!data) {
        fprintf(stderr, "Cannot load ImageData \"%s\"\nSTB Reason: %s\n", filename, stbi_failure_reason());
        exit(0);
    }
    if (channels) c = channels;
    int i, j, k;
    image im = make_image(h, w, c);
    for (k = 0; k < c; ++k) {
        for (j = 0; j < h; ++j) {
            for (i = 0; i < w; ++i) {
                int dst_index = i + w * j + w * h*k;
                int src_index = k + c * i + c * w*j;
                im.data[dst_index] = (float)data[src_index] / 255.;
            }
        }
    }
    free(data);
    return im;
}

image load_image(char *filename)
{
    image out = load_image_stb(filename, 3);
    return out;
	
	/*
    IplImage* src = 0;
    if( (src = cvLoadImage(filename,-1)) == 0 )
    {
        printf("Cannot load file image %s\n", filename);
        exit(0);
    }
    unsigned char *data = (unsigned char *)src->imageData;
    int c = src->nChannels;
    int h = src->height;
    int w = src->width;
    int step = src->widthStep;
    image out = make_image(h,w,c);
    int i, j, k, count=0;;

    for(k= 0; k < c; ++k){
        for(i = 0; i < h; ++i){
            for(j = 0; j < w; ++j){
                out.data[count++] = data[i*step + j*c + k];
            }
        }
    }
    cvReleaseImage(&src);
    return out;
	*/
}

image get_image_layer(image m, int l)
{
    image out = make_image(m.h, m.w, 1);
    int i;
    for(i = 0; i < m.h*m.w; ++i){
        out.data[i] = m.data[i+l*m.h*m.w];
    }
    return out;
}

double get_pixel(image m, int x, int y, int c)
{
    assert(x < m.h && y < m.w && c < m.c);
    return m.data[c*m.h*m.w + x*m.w + y];
}
double get_pixel_extend(image m, int x, int y, int c)
{
    if(x < 0 || x >= m.h || y < 0 || y >= m.w || c < 0 || c >= m.c) return 0;
    return get_pixel(m, x, y, c);
}
void set_pixel(image m, int x, int y, int c, double val)
{
    assert(x < m.h && y < m.w && c < m.c);
    m.data[c*m.h*m.w + x*m.w + y] = val;
}
void set_pixel_extend(image m, int x, int y, int c, double val)
{
    if(x < 0 || x >= m.h || y < 0 || y >= m.w || c < 0 || c >= m.c) return;
    set_pixel(m, x, y, c, val);
}

void add_pixel(image m, int x, int y, int c, double val)
{
    assert(x < m.h && y < m.w && c < m.c);
    m.data[c*m.h*m.w + x*m.w + y] += val;
}

void add_pixel_extend(image m, int x, int y, int c, double val)
{
    if(x < 0 || x >= m.h || y < 0 || y >= m.w || c < 0 || c >= m.c) return;
    add_pixel(m, x, y, c, val);
}

void two_d_convolve(image m, int mc, image kernel, int kc, int stride, image out, int oc, int edge)
{
    int x,y,i,j;
    int xstart, xend, ystart, yend;
    if(edge){
        xstart = ystart = 0;
        xend = m.h;
        yend = m.w;
    }else{
        xstart = kernel.h/2;
        ystart = kernel.w/2;
        xend = m.h-kernel.h/2;
        yend = m.w - kernel.w/2;
    }

    for(x = xstart; x < xend; x += stride)	{
        for(y = ystart; y < yend; y += stride){
            double sum = 0;
            for(i = 0; i < kernel.h; ++i){
                for(j = 0; j < kernel.w; ++j){
                    sum += get_pixel(kernel, i, j, kc)*get_pixel_extend(m, x+i-kernel.h/2, y+j-kernel.w/2, mc);
                }
            }
            add_pixel(out, (x-xstart)/stride, (y-ystart)/stride, oc, sum);
        }
    }
}

double single_convolve(image m, image kernel, int x, int y)
{
    double sum = 0;
    int i, j, k;
    for(i = 0; i < kernel.h; ++i){
        for(j = 0; j < kernel.w; ++j){
            for(k = 0; k < kernel.c; ++k){
                sum += get_pixel(kernel, i, j, k)*get_pixel_extend(m, x+i-kernel.h/2, y+j-kernel.w/2, k);
            }
        }
    }
    return sum;
}

void convolve(image m, image kernel, int stride, int channel, image out, int edge)
{
    assert(m.c == kernel.c);
    int i;
    zero_channel(out, channel);
    for(i = 0; i < m.c; ++i)
	{
        two_d_convolve(m, i, kernel, i, stride, out, channel, edge);
    }
    /*
    int j;
    for(i = 0; i < m.h; i += stride){
        for(j = 0; j < m.w; j += stride){
            double val = single_convolve(m, kernel, i, j);
            set_pixel(out, i/stride, j/stride, channel, val);
        }
    }
    */
}

void upsample_image(image m, int stride, image out)
{
    int i,j,k;
    zero_image(out);
    for(k = 0; k < m.c; ++k){
        for(i = 0; i < m.h; ++i){
            for(j = 0; j< m.w; ++j){
                double val = get_pixel(m, i, j, k);
                set_pixel(out, i*stride, j*stride, k, val);
            }
        }
    }
}

void single_update(image m, image update, int x, int y, double error)
{
    int i, j, k;
    for(i = 0; i < update.h; ++i){
        for(j = 0; j < update.w; ++j){
            for(k = 0; k < update.c; ++k){
                double val = get_pixel_extend(m, x+i-update.h/2, y+j-update.w/2, k);
                add_pixel(update, i, j, k, val*error);
            }
        }
    }
}

void kernel_update(image m, image update, int stride, int channel, image out, int edge)
{
    assert(m.c == update.c);
    zero_image(update);
    int i, j, istart, jstart, iend, jend;
    if(edge){
        istart = jstart = 0;
        iend = m.h;
        jend = m.w;
    }else{
        istart = update.h/2;
        jstart = update.w/2;
        iend = m.h-update.h/2;
        jend = m.w - update.w/2;
    }
    for(i = istart; i < iend; i += stride){
        for(j = jstart; j < jend; j += stride){
            double error = get_pixel(out, (i-istart)/stride, (j-jstart)/stride, channel);
            single_update(m, update, i, j, error);
        }
    }
    /*
    for(i = 0; i < update.h*update.w*update.c; ++i){
        update.data[i] /= (m.h/stride)*(m.w/stride);
    }
    */
}

void single_back_convolve(image m, image kernel, int x, int y, double val)
{
    int i, j, k;
    for(i = 0; i < kernel.h; ++i){
        for(j = 0; j < kernel.w; ++j){
            for(k = 0; k < kernel.c; ++k){
                double pval = get_pixel(kernel, i, j, k) * val;
                add_pixel_extend(m, x+i-kernel.h/2, y+j-kernel.w/2, k, pval);
            }
        }
    }
}

void back_convolve(image m, image kernel, int stride, int channel, image out, int edge)
{
    assert(m.c == kernel.c);
    int i, j, istart, jstart, iend, jend;
    if(edge){
        istart = jstart = 0;
        iend = m.h;
        jend = m.w;
    }else{
        istart = kernel.h/2;
        jstart = kernel.w/2;
        iend = m.h-kernel.h/2;
        jend = m.w - kernel.w/2;
    }
    for(i = istart; i < iend; i += stride){
        for(j = jstart; j < jend; j += stride){
            double val = get_pixel(out, (i-istart)/stride, (j-jstart)/stride, channel);
            single_back_convolve(m, kernel, i, j, val);
        }
    }
}

void print_image(image m)
{
    int i;
    for(i =0 ; i < m.h*m.w*m.c; ++i) printf("%lf, ", m.data[i]);
    printf("\n");
}

void free_image(image m)
{
    free(m.data);
}
