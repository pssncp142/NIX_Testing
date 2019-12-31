#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "cpl_image.h"
#include "cpl.h"
#include "hdrl.h"

bool ERROR_IMAGE = false;
int ERROR_METHOD = 0;
double GAIN = 5.7;
double RN_ADU = 4.5;


void hdrl_init(){
	cpl_init(CPL_INIT_DEFAULT);
}

void hdrl_end(){
	cpl_end();
}

hdrl_parameter * hdrl_bpm_fit_parameter_create(int degree, double pval, double rel_chi_l, double rel_chi_h, double rel_coef_l, double rel_coef_h) {

	if (pval != -1){
		return hdrl_bpm_fit_parameter_create_pval(degree, pval);
	} else if (rel_chi_l != -1) {
		return hdrl_bpm_fit_parameter_create_rel_chi(degree, rel_chi_l, rel_chi_h);
	} else if (rel_coef_l != -1) {
		return hdrl_bpm_fit_parameter_create_rel_coef(degree, rel_coef_l, rel_coef_h);
	}

}

hdrl_parameter * hdrl_bpm_2d_parameter_create(cpl_filter_mode filter, cpl_border_mode border,
											double kappa_low, double kappa_high, int maxiter,
											int steps_x, int steps_y,
											int filter_size_x, int filter_size_y,
											int order_x, int order_y,
											int smooth_x, int smooth_y){

	printf("%d\n", order_x);
	if (order_x == -1){
		return hdrl_bpm_2d_parameter_create_filtersmooth(kappa_low, kappa_high, maxiter, filter, border, smooth_x, smooth_y);
	} else {
		return hdrl_bpm_2d_parameter_create_legendresmooth(kappa_low, kappa_high, maxiter,
															steps_x, steps_y,
															filter_size_x, filter_size_y,
															order_x, order_y);

	}

}

hdrl_image * hdrl_image_create_numpy_float64(double * image, int nx, int ny){

	cpl_image * image_cpl, * image_err;

	image_cpl = cpl_image_wrap_double(nx, ny, image);

	if (ERROR_IMAGE == true) {

		if (ERROR_METHOD == 0){
			image_err = cpl_image_abs_create(image_cpl);
			cpl_image_multiply_scalar(image_err, GAIN);
			cpl_image_add_scalar(image_err, RN_ADU*RN_ADU*GAIN*GAIN);
			cpl_image_power(image_err, 0.5);
			cpl_image_divide_scalar(image_err, GAIN);
		} else {
			image_err = cpl_image_new(nx, ny, CPL_TYPE_DOUBLE);
			cpl_image_add_scalar(image_err, cpl_image_get_median(image_cpl)*GAIN);
			cpl_image_power(image_err, 0.5);
			cpl_image_divide_scalar(image_err, GAIN);
		}
	} else {
		image_err = cpl_image_new(nx, ny, CPL_TYPE_DOUBLE);
	}

	hdrl_image * image_hdrl;

	image_hdrl = hdrl_image_create(image_cpl, image_err);

	cpl_image_unwrap(image_cpl);
	cpl_image_unwrap(image_err);

	return image_hdrl;

}

hdrl_imagelist * hdrl_imagelist_create_numpy_float64(double * image, int nx, int ny, int nz){

	hdrl_imagelist * image_list;
	hdrl_image * image_hdrl;

	image_list = hdrl_imagelist_new();

	for (int i=0; i<nz; i++){
		image_hdrl = hdrl_image_create_numpy_float64(&image[nx*ny*i], nx, ny);
		hdrl_imagelist_set(image_list, image_hdrl, i);
	}


	return image_list;

}

cpl_vector * cpl_vector_create_numpy_float64(double * items, int n){

	cpl_vector * vector;

	vector = cpl_vector_new(n);

	for (int i=0; i<n; i++){
		cpl_vector_set(vector, i, items[i]);
	}

	return vector;

}

hdrl_strehl_result hdrl_compute_strehl_numpy_float64(double * image, int nx, int ny,
		double wave, double m1_rad, double m2_rad,
		double pixsc_x, double pixsc_y,
		double flux_r, double bkg_low, double bkg_high){

	hdrl_image * image_hdrl;
	hdrl_strehl_result strehl;

	image_hdrl = hdrl_image_create_numpy_float64(image, nx, ny);

	hdrl_parameter * params = hdrl_strehl_parameter_create(wave, m1_rad, m2_rad,
							pixsc_x, pixsc_y,
							flux_r, bkg_low, bkg_high);

	strehl = hdrl_strehl_compute(image_hdrl, params);

	hdrl_image_delete(image_hdrl);
	hdrl_parameter_delete(params);

	return strehl;

}

int hdrl_bpm_2d_compute_numpy_float64(double * image, cpl_binary * mask_out, hdrl_parameter * params, int nx, int ny){

	hdrl_image * image_hdrl;
	cpl_mask * bp_mask;
	int result;

	image_hdrl = hdrl_image_create_numpy_float64(image, nx, ny);
	bp_mask =  hdrl_bpm_2d_compute(image_hdrl, params);
	memcpy(mask_out, cpl_mask_get_data(bp_mask), sizeof(cpl_binary)*nx*ny);

	hdrl_image_delete(image_hdrl);

	return result;

}

int hdrl_bpm_3d_compute_numpy_float64(double * images_in, cpl_binary * masks_out, int nx, int ny, int nz){
	
	hdrl_imagelist * image_list;
	cpl_imagelist * bpm_3d;
	cpl_image * image_cpl;
	int * tmp_mask;

	image_list = hdrl_imagelist_create_numpy_float64(images_in, nx, ny, nz);

	//HDRL_BPM_3D_THRESHOLD_ABSOLUTE
	//HDRL_BPM_3D_THRESHOLD_RELATIVE
	//HDRL_BPM_3D_THRESHOLD_ERROR

	hdrl_parameter * params = hdrl_bpm_3d_parameter_create(5, 10, HDRL_BPM_3D_THRESHOLD_RELATIVE);
	bpm_3d =  hdrl_bpm_3d_compute(image_list, params);

	for (int i=0; i<nz; i++){
		image_cpl = cpl_imagelist_get(bpm_3d, i);
		tmp_mask = (int*) cpl_image_get_data(image_cpl);
		for(int j=0; j<nx*ny; j++){
			masks_out[i*nx*ny+j] = (cpl_binary) tmp_mask[j];
		}
	}

	hdrl_imagelist_delete(image_list);
	hdrl_parameter_delete(params);
	cpl_imagelist_delete(bpm_3d);

	return 0;

}

int hdrl_bpm_fit_compute_numpy_float64(double * images_in, cpl_binary * mask_out, double * exptime, hdrl_parameter * params, int nx, int ny, int nz){
	
	//cpl_init(CPL_INIT_DEFAULT);

	hdrl_imagelist * image_list;
	cpl_vector *  sample_pos;
	cpl_image * mask;
	int * tmp_mask;

	image_list = hdrl_imagelist_create_numpy_float64(images_in, nx, ny, nz);

	sample_pos = cpl_vector_create_numpy_float64(exptime, nz);

	int res = hdrl_bpm_fit_compute(params, image_list, sample_pos, &mask);

	tmp_mask = (int*) cpl_image_get_data(mask);

	for (int i=0; i<nx*ny; i++){
		mask_out[i] = (cpl_binary) tmp_mask[i];
	}

	
	hdrl_imagelist_delete(image_list);
	//hdrl_parameter_delete(params);

	//cpl_end();

	return 0;

}

int hdrl_collapse_median_numpy_float64(double * images_in, double * image_out, int nx, int ny, int nz){

	hdrl_imagelist * image_list;
	hdrl_image * image_hdrl;
	cpl_image * image_cpl;
	hdrl_image * result; cpl_image * contrib;
	double * image_p;

	image_list = hdrl_imagelist_create_numpy_float64(images_in, nx, ny, nz);

	int error = hdrl_imagelist_collapse(image_list, HDRL_COLLAPSE_MEDIAN, &result, &contrib); 
	image_cpl = hdrl_image_get_image(result);
	image_p = (double*) cpl_image_get_data(image_cpl);
	memcpy(image_out, image_p, sizeof(double)*nx*ny);

	hdrl_imagelist_delete(image_list);
	cpl_image_unwrap(image_cpl);
	free(image_p);
	
	return 0;

}


