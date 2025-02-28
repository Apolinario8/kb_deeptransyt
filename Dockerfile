FROM kbase/sdkpython:3.8.0
MAINTAINER goncalo_apolinario
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# RUN apt-get update


# -----------------------------------------
# RUN apt-get update && \
#     apt-get install -y wget bzip2 && \
#     wget -nv https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
#     bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
#     rm Miniconda3-latest-Linux-x86_64.sh

# ENV PATH=/opt/conda/bin:$PATH

# RUN conda create -n deeptransyt_env python=3.8

# RUN /opt/conda/bin/conda install -n deeptransyt_env pip && \
#     /opt/conda/envs/deeptransyt_env/bin/pip install --upgrade pip && \
#     /opt/conda/envs/deeptransyt_env/bin/pip install deeptransyt==0.0.3

# COPY ./ /kb/module
# RUN mkdir -p /kb/module/work
# RUN chmod -R a+rw /kb/module

# WORKDIR /kb/module

# RUN make all

# #ENTRYPOINT [ "./scripts/entrypoint.sh" ]
# ENTRYPOINT [ "bash", "-c", "source activate deeptransyt_env && ./scripts/entrypoint.sh" ]

# CMD [ ]

RUN apt-get update && \
    apt-get install -y wget bzip2 && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip
RUN pip install deeptransyt==0.0.13  
COPY ./ /kb/module

RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "bash", "-c", "./scripts/entrypoint.sh" ]

CMD [ ]
