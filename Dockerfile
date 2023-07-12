FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN pip install rdkit==2023.3.2 
RUN pip install mordred==1.2
RUN pip install timeout-decorator==0.5.0
RUN pip install pandas==1.3.5

WORKDIR /repo
COPY ./repo