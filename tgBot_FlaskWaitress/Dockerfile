FROM python:3.7.9-stretch

RUN apt-get install \
        ca-certificates \
        apt-transport-https \
        curl

RUN pip install pipenv

WORKDIR /app

COPY . /app

ENTRYPOINT ["pipenv", "run", "python", "configure_and_run.py", "&"]
