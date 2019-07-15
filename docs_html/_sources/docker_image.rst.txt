Docker image
============

`Docker <https://docker.com/>`_ is a light-weight way to deploy
applications. We offer an `image for READemption
<https://registry.hub.docker.com/u/konradfoerstner/reademption/>`_
that can be easily retrieved and used once you have `Docker installed
<https://docs.docker.com/installation/>`_. To get the image run

::

  $ sudo docker pull konradfoerstner/reademption

Now you have the image READempmtion image and should see it shown
in your list of docker images:

::

  $ sudo docker images
  REPOSITORY                    TAG                 IMAGE ID            CREATED             VIRTUAL SIZE
  reademption_0.3.3             latest              416dbb001c6b        2 weeks ago         802.5 MB
  konradfoerstner/reademption   0.3.3               416dbb001c6b        2 weeks ago         802.5 MB
  ubuntu                        14.04               ba5877dc9bec        4 weeks ago         192.7 MB
  ubuntu                        latest              ba5877dc9bec        4 weeks ago         192.7 MB

You can simply start the contains by running the following command:

:: 

  $ sudo docker run -i -t reademption_0.3.3 
  root@3b9e7fd860ee:/# reademption -v
  READemption version 0.3.3

Inside of this environment you have access to the local installation
of READemption.
