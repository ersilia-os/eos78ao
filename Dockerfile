FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN pip install rdkit==2023.3.2 
RUN pip install mordred==1.2
RUN pip install timeout-decorator==0.5.0
RUN pip install pandas==1.3.5
RUN conda install -c conda-forge xorg-libxrender xorg-libxtst

WORKDIR /repo
COPY ./repo
