#Dockerfile for our program
FROM python:3.7
ADD src/primerPipe.py .
RUN pip install -r requirements.txt
