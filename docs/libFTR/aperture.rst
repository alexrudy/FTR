.. highlight:: c
.. default-domain:: c

Aperture management tools
*************************

There are few simple aperture managment tools in ``aperture.h`` which might be useful for some simple applications. In general, it is expected that the user specify the aperture very precisely in an integer array, with ``1`` for illuminated subapertures and ``0`` for non-illuminated subapertures. This set of routines is generally for circular apertures, though the user could provide their own aperture array and overwrite the interal aperture arry setting.

Reference / API
===============

The :type:`aperture` structure
------------------------------

.. type:: aperture

    This is the basic aperture type. It is a pointer to the aperture structure type.

.. type:: aperture_s

    The raw aperture structure.

.. member:: aperture_s.nx

    Number of subapertures in the x direction.

.. member:: aperture_s.ny

    Number of subapertures in the y direction.

.. member:: aperture_s.nm

    Number of modes controlled for this aperture. This value is normally set to 0 by code in this header, but can be useful for other user-facing features.

.. member:: aperture_s.ap

    A pointer to the aperture illumination array, with ``1`` for illuminated subapertures and ``0`` for non-illuminated subapertures.

.. member:: aperture_s.ni

    The number of illuminated subapertures.

Creating apertures
------------------

.. function:: aperture aperture_create(const int ny, const int nx, const int *ap)

    Create an aperture from an integer array of illiuminated points.
    
    :param const int ny: Number of points in the Y direction.
    :param const int nx: Number of points in the X direction.
    :param const int* ap: Aperture array.
    :returns :type:`aperture` ap: A pointer to an aperture structure.


.. function:: void aperture_destroy(aperture ap)
    
    Destroy and deallocate aperture arrays.
    
    :param :type:`aperture` ap: Aperture object to deallocate.
    

.. function:: aperture aperture_create_default(const int ny, const int nx)

    This creates an aperture with a default annulus pattern that almost fills the entire map, and has a secondary obscuration.
    It is good for rough testing. The outer radius is set by ``(n / 2) - 1``, and the inner radius is set to ``outer_radius / 3``.
    
    :param const int ny: Number of points in the Y direction.
    :param const int nx: Number of points in the X direction.
    :returns :type:`aperture` ap: A pointer to an aperture structure.
    

.. function:: aperture aperture_create_with_radii(const int ny, const int nx, const double outer_radius, const double inner_radius)
    
    This function creates an aperture from a cicular annulus with an inner and outer radius.
    
    :param const int ny: Number of points in the Y direction.
    :param const int nx: Number of points in the X direction.
    :param const double outer_radius: The outer radius of the annulus.
    :param const double inner_radius: The inner radius of the annulus.
    :returns :type:`aperture` ap: A pointer to an aperture structure.

Viewing apertures
-----------------

.. function:: void aperture_print(const aperture ap)
    
    Print the aperture mask to stdout.
    
    :param :type:`aperture` ap: Aperture to print.
    
    


