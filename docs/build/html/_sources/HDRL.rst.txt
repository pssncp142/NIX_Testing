NIX_Testing.HDRL2
=================

.. toctree::
    :maxdepth: 2

To simplify, CPL images are created of type double and cpl masks of type cpl_binary (or unsigned char).
Their equivalents in numpy are NPY_FLOAT64 and NPY_UBYTE respectively.
The functions in this module accept ndarrays of any type, but only returns the types given above.
Any type beyond np.float64 will potentially throw an error (due to explicit casting safety on numpy), or potentially segfault if given as input.
Therefore, they should be avoided!.

Global Variables
----------------

HDRL image struct is created with data and its error. Some functions in the library make use of the error values. 
Hence, this module generates error images from the data based on the global variables given below. 

.. data:: NIX_Testing.HDRL2.GAIN=5.7 (*float*)
.. data:: NIX_Testing.HDRL2.RN_ADU=4.5 (*float*)
.. data:: NIX_Testing.HDRL2.ERROR_IMAGE=True (*bool*)
.. data:: NIX_Testing.HDRL2.ERROR_METHOD=0 (*int*)

Generate Bad Pixel Maps
-----------------------

.. autofunction:: NIX_Testing.HDRL2.bpm_fit_compute

.. automethod:: NIX_Testing.HDRL2.bpm_2d_compute

Cosmic-rays
-----------

.. automethod:: NIX_Testing.HDRL2.lacosmic_edgedetect

Interpolate
-----------

.. automethod:: NIX_Testing.HDRL2.bpm_interpolate

Image Quality
-------------

.. automethod:: NIX_Testing.HDRL2.compute_strehl

