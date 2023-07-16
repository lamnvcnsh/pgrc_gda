# Use an official Ubuntu runtime as a parent image
FROM ubuntu:20.04

# Avoid timezone interactive dialog during installation
ENV DEBIAN_FRONTEND=noninteractive

# Install Python and pip
RUN apt-get update && apt-get install -y python3.8 python3-pip openjdk-11-jdk
RUN apt-get install -y libbz2-dev liblzma-dev

# Alias python3.8 to python
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.8 1

# Set the working directory in the container to /app
WORKDIR /app

# Add the current directory contents into the container at /app
ADD . /app

RUN mkdir data
RUN chmod 755 run_stargazer.sh data stargazer

# Upgrade pip and install Python dependencies
RUN pip install pysam pandas scipy matplotlib

# Run your application when the container launches
# CMD [ "python setup.py install" ]
