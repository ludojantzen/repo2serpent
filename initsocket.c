/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : initsocket.c                                   */
/*                                                                           */
/* Created:       2018/03/15 (VVa)                                           */
/* Last modified: 2020/05/28 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Initializes the socket and communicates to the server.       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InitSocket:"

/* MACOS compability */

#define h_addr h_addr_list[0]

/*****************************************************************************/

void InitSocket()
{
  long port, my_socket, i;
  struct sockaddr_in serv_addr;
  struct hostent *server;
  char buffer[MAX_STR];

  /* Return if no parent port has been given */

  if ((port = (long)RDB[DATA_PPORT]) == 0)
    return;

  /* Open a socket (IP-stream) */

  my_socket = socket(AF_INET, SOCK_STREAM, 0);

  /* Check that the socket was opened */

  if (my_socket < 0)
    Die(FUNCTION_NAME, "Could not open socket.");

  /* Currently our server will always be on the same host as Serpent */

  server = gethostbyname("localhost");

  /* This should not fail */

  if (server == NULL)
    Die(FUNCTION_NAME,"ERROR, no such host\n");

  /* Zero out the server address structure */

  bzero((char *) &serv_addr, sizeof(serv_addr));

  /* Fill in the known information about the server address */

  serv_addr.sin_family = AF_INET;
  bcopy((char *)server->h_addr, (char *)&serv_addr.sin_addr.s_addr,
        server->h_length);

  for (i = 19; i >= 0; i--)
    {
      /* Try next port */

      serv_addr.sin_port = htons(port+i);

      /* Connect to the server */

      if (connect(my_socket, (struct sockaddr *)&serv_addr,sizeof(serv_addr)) < 0)
        fprintf(outp, "Could not connect to port %ld, trying the next one.\n", port+i);
      else
        break;
    }

  if (i == 20)
    Die(FUNCTION_NAME, "ERROR connecting to ports %ld--%ld", port, port+19);

  /* Store socket */

  WDB[DATA_COM_SOCKET] = (double)my_socket;

  /* Send a test string */

  sprintf(buffer, "Serpent has connected");
  SocketSendString(buffer);

  /* Zero out the buffer */

  bzero(buffer, MAX_STR);

  /* Receive a test string */

  SocketReceiveString(buffer);

  /* Print out the test string */

  fprintf(outp, "%s\n",buffer);

  return;

}

/*****************************************************************************/
