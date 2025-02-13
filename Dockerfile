FROM bentoml/model-server:0.11.0-py312
MAINTAINER ersilia

RUN pip install mordredcommunity[full]==2.0.6
RUN pip install timeout-decorator==0.5.0
RUN pip install scikit-learn==1.6.1
RUN pip install joblib==1.4.2

WORKDIR /repo
COPY ./repo
