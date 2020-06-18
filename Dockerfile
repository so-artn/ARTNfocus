FROM python:latest

MAINTAINER T. E. Pickering "te.pickering@gmail.com"

COPY . .

RUN pip install -e .

WORKDIR /artn

ENTRYPOINT ["mont4k_focus"]
