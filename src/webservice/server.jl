server = nothing
port = 8000

function start_server(_port::Int64)
	log_message("starting up server...")
	global port = _port
	global server = Sockets.listen(Sockets.InetAddr(Sockets.localhost, port))
    # global server = HTTP.Servers.Server(create_base_router())
    try
        HTTP.serve(create_base_router(), "127.0.0.1", port, verbose = true, server = server)
    finally
        close(server)
	end
end

stop_server() = close(server)
