Python module
===================================
.. role:: python(code)
   :language: python

.. py:module:: fstpack
.. py:function:: fst (a)

   Compute the one-dimensional discrete Stockwell Transform.

   This function computes the one-dimensional discrete Stockwell Transform with the efficient Fast Stockwell Transform (FST) algorithm [ST]_.

   :param array_like a: Input array, can be complex

   :return: The transformed input
   :rtype: complex ndarray

.. py:function:: ifst(a)

   Compute the one-dimensional inverse discrete Stockwell Transform.

   This function computes the inverse of the one-dimensional discrete Stockwell transform computed by :py:func:`fst`.
   In other words, :python:`ifst(fst(a)) == a` to within numerical accuracy.

   :param array_like a: Input array, can be complex

   :return: The transformed input
   :rtype: complex ndarray

.. py:function:: dst2(a)

.. py:function:: idst2(a)

.. py:function:: imfreq(a, x, y)

Notes
-------------

References
-------------
.. [ST] Article

Examples
-------------
.. code-block:: python

   >>> np.fft.fft(np.exp(2j * np.pi * np.arange(8) / 8))
   array([-2.33486982e-16+1.14423775e-17j,  8.00000000e+00-1.25557246e-15j,
        2.33486982e-16+2.33486982e-16j,  0.00000000e+00+1.22464680e-16j,
       -1.14423775e-17+2.33486982e-16j,  0.00000000e+00+5.20784380e-16j,
        1.14423775e-17+1.14423775e-17j,  0.00000000e+00+1.22464680e-16j])
