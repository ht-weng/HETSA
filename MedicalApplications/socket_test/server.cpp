#include <stdio.h>
#include <unistd.h>
#include <netinet/in.h>
#include <string.h>
#include <sys/socket.h>
#include <stdlib.h>
#include <iostream>
#define PORT 8080
using namespace std;

int main(int argument, char const *argv[])
{
    int obj_server, sock, reader;
    struct sockaddr_in address;
    int opted = 1;
    int address_length = sizeof(address);
    char buffer[1024] = {0};
    const char *message = "A message from server !";
    if ((obj_server = socket(AF_INET, SOCK_STREAM, 0)) == 0)
    {
        cout << "Opening of Socket Failed !" << endl;
        exit(EXIT_FAILURE);
    }
    if (setsockopt(obj_server, SOL_SOCKET, SO_REUSEADDR,
                     &opted, sizeof(opted)))
    {
        cout << "Can't set the socket" << endl;
        exit(EXIT_FAILURE);
    }
    address.sin_family = AF_INET;
    address.sin_addr.s_addr = INADDR_ANY;
    address.sin_port = htons(PORT);
    if (bind(obj_server, (struct sockaddr *)&address,
             sizeof(address)) < 0)
    {
        cout << "Binding of socket failed !" << endl;
        exit(EXIT_FAILURE);
    }
    if (listen(obj_server, 3) < 0)
    {
        cout << "Can't listen from the server !" << endl;
        exit(EXIT_FAILURE);
    }
    if ((sock = accept(obj_server, (struct sockaddr *)&address, (socklen_t *)&address_length)) < 0)
    {
        cout << "Accept" << endl;
        exit(EXIT_FAILURE);
    }
    reader = read(sock, buffer, 1024);
    cout << buffer << endl;
    send(sock, message, strlen(message), 0);
    cout << "Server : Message has been sent ! \n" << endl;
    return 0;
}