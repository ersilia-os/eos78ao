FROM bentoml/model-server:0.11.0-py38
MAINTAINER ersilia

RUN pip install mordredcommunity==2.0.6
RUN pip install networkx==3.2.1
RUN pip install numpy==1.26.4
RUN pip install pandas==1.3.5
RUN pip install rdkit==2023.3.2
RUN pip install timeout-decorator==0.5.0

WORKDIR /repo
COPY ./repo
