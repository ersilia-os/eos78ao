FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN conda install -c conda-forge rdkit=2021.03
RUN conda install -c mordred-descriptor mordred=1.2
RUN pip install timeout-decorator==0.5.0

WORKDIR /repo
COPY ./repo
