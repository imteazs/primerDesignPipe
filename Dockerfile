#Dockerfile for our program
FROM python:3.7.12-bullseye
WORKDIR /primerdesign
COPY requirements.txt .
RUN pip install -r requirements.txt
COPY src/primerPipe.py .
CMD ["python", "primerPipe.py", "-h"]
