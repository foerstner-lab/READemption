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
  tillsauerwein/reademption     2.0.2               41c42c0a6c35        42 minutes ago      3.31GB
  tillsauerwein/reademption     latest              241c42c0a6c35        42 minutes ago      3.31GB



You can simply start the container by running the following command:

:: 

  $ sudo docker run -i -t reademption
  root@3b9e7fd860ee:/# reademption -v
  READemption version 2.0.2


Inside of this environment you have access to the local installation
of READemption.

You can find a detailed tutorial of using READemption with docker at https://github.com/Tillsa/READemption_Docker_Tutorial.
The tutorial is the same as the single-species analysis shown
at the `Performing example analyses <https://reademption.readthedocs.io/en/latest/example_analysis.html#single-species-analysis>`_ section.
