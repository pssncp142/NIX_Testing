#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cpl_image.h"
#include "cpl.h"
#include "hdrl.h"



hdrl_image * hdrl_image_create_numpy_float64(double * image, int nx, int ny){

	cpl_image * image_cpl, * image_err;

	image_cpl = cpl_image_wrap_double(nx, ny, image);
	image_err = cpl_image_new(nx, ny, CPL_TYPE_DOUBLE);

	cpl_image_add_scalar(image_err, cpl_image_get_median(image_cpl)/5.);
	cpl_image_power(image_err, 0.5);

	//cpl_image_multiply_scalar(image_cpl, 5.);
	//cpl_image_copy(image_err, image_cpl, 1, 1);
	//cpl_image_power(image_err, 0.5);

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

	cpl_init(CPL_INIT_DEFAULT);

	hdrl_image * image_hdrl;
	hdrl_strehl_result strehl;

	image_hdrl = hdrl_image_create_numpy_float64(image, nx, ny);

	hdrl_parameter * params = hdrl_strehl_parameter_create(wave, m1_rad, m2_rad,
							pixsc_x, pixsc_y,
							flux_r, bkg_low, bkg_high);

	strehl = hdrl_strehl_compute(image_hdrl, params);

	hdrl_image_delete(image_hdrl);
	hdrl_parameter_delete(params);

	cpl_end();

	return strehl;

}

int hdrl_bpm_2d_compute_numpy_float64(double * image, cpl_binary * mask_out, int nx, int ny){

	cpl_init(CPL_INIT_DEFAULT);

	hdrl_image * image_hdrl;
	cpl_mask * bp_mask;
	int result;

	image_hdrl = hdrl_image_create_numpy_float64(image, nx, ny);

	hdrl_parameter * params =  hdrl_bpm_2d_parameter_create_filtersmooth(5., 10., 10, CPL_FILTER_MEDIAN, CPL_BORDER_FILTER, 5, 5);
	//hdrl_parameter * params =  hdrl_bpm_2d_parameter_create_legendresmooth(5., 10., 10, 10, 10, 10, 10, 2, 2);

	bp_mask =  hdrl_bpm_2d_compute(image_hdrl, params);

	result = cpl_mask_save(bp_mask, f_name, NULL, CPL_IO_CREATE);

	memcpy(mask_out, cpl_mask_get_data(bp_mask), sizeof(cpl_binary)*nx*ny);

	hdrl_image_delete(image_hdrl);
	hdrl_parameter_delete(params);

	cpl_end();

	return result;

}

int hdrl_bpm_3d_compute_numpy_float64(double * images_in, cpl_binary * masks_out, int nx, int ny, int nz){
	
	cpl_init(CPL_INIT_DEFAULT);

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

	cpl_end();

	return 0;

}

int hdrl_bpm_fit_compute_numpy_float64(double * images_in, cpl_binary * mask_out, double * exptime, int nx, int ny, int nz){
	
	cpl_init(CPL_INIT_DEFAULT);

	hdrl_imagelist * image_list;
	cpl_vector *  sample_pos;
	cpl_image * mask;
	int * tmp_mask;

	image_list = hdrl_imagelist_create_numpy_float64(images_in, nx, ny, nz);

	hdrl_parameter * params = hdrl_bpm_fit_parameter_create_pval(3., 20.);

	sample_pos = cpl_vector_create_numpy_float64(exptime, nz);

	int res = hdrl_bpm_fit_compute(params, image_list, sample_pos, &mask);

	tmp_mask = (int*) cpl_image_get_data(mask);

	for (int i=0; i<nx*ny; i++){
		mask_out[i] = (cpl_binary) tmp_mask[i];
	}

	
	hdrl_imagelist_delete(image_list);
	hdrl_parameter_delete(params);

	cpl_end();

	return 0;

}

int hdrl_collapse_median_numpy_float64(double * images_in, double * image_out, int nx, int ny, int nz){

	cpl_init(CPL_INIT_DEFAULT);

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

	cpl_end();
	
	return 0;

}


