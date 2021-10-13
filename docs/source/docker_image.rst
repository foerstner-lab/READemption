Docker image
============

`Docker <https://docker.com/>`_ is a light-weight way to deploy
applications. We offer an `image for READemption
<https://hub.docker.com/r/tillsauerwein/reademption>`_
that can be easily retrieved and used once you have `Docker installed
<https://docs.docker.com/installation/>`_. To get the image run

::

  $ sudo docker pull tillsauerwein/reademption

Now you have the image READempmtion image and should see it shown
in your list of docker images:

::

  $ sudo docker images
  REPOSITORY                    TAG                 IMAGE ID            CREATED             VIRTUAL SIZE
  tillsauerwein/reademption     1.0.5               6ca6b40723c0        3 hours ago         1.91GB
  tillsauerwein/reademption     latest              6ca6b40723c0        3 hours ago         1.91GB


You can simply start the contains by running the following command:

:: 

  $ sudo docker run -i -t reademption_1.0.5
  root@3b9e7fd860ee:/# reademption -v
  READemption version 1.0.5


Inside of this environment you have access to the local installation
of READemption.

You can find a detailed tutorial of using READemption with docker at https://github.com/Tillsa/READemption_Docker_Tutorial.
