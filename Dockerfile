#Dockerfile for our program
FROM python:3.7.12-bullseye
ADD src/primerPipe.py .
RUN pip install -r requirements.txt
