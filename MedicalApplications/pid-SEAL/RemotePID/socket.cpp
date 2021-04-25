#include <iostream>
#include <ctime>
#include <string>
// #include <boost/bind.hpp>
// #include <boost/shared_ptr.hpp>
// #include <boost/enable_shared_from_this.hpp>
#include <boost/asio.hpp>

using namespace boost::asio;
using namespace std;

int main() {
    ip::tcp::iostream stream;

    if (!stream) {
        cout << "Can't connect" << endl;
    }

    // stream.connect("127.0.0.1", "http");
    // stream << "GET /LICENSE_1_0.txt HTTP/1.0\r\n";
    // stream << "Host: www.boost.org\r\n";
    // stream << "Accept: */*\r\n";
    // stream << "Connection: close\r\n\r\n";
    // stream.flush();
    // cout << stream.rdbuf();

    return 0;
}

