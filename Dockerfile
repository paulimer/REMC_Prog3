FROM registry.fedoraproject.org/fedora:36
MAINTAINER Paul Etheimer <paul.etheimer@etu.u-paris.fr>
RUN dnf -y upgrade && dnf install -y \
    python3-numpy \
    && dnf clean all
RUN groupadd --system remc && adduser --system --gid remc remc
WORKDIR /home/remc
COPY . .
RUN chown remc:remc -R /home/remc
USER remc
RUN chmod u+x ./bin/main.py
