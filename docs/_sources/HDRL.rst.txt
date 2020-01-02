NIX_Testing.HDRL2
=================

.. toctree::
    :maxdepth: 2

.. automodule:: NIX_Testing.HDRL2
    :members:

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

