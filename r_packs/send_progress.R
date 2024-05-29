library(reticulate)

# 定义Python代码，用于发送进度信息到RabbitMQ
py_run_string("
import pika

def send_progress(message, queue='main_queue'):
    connection = pika.BlockingConnection(pika.ConnectionParameters('localhost'))
    channel = connection.channel()
    channel.queue_declare(queue)
    channel.basic_publish(exchange='', routing_key=queue, body=message)
    connection.close()
")

send_progress <- import_main()$send_progress

progress_callback <- function(message) {
    send_progress(message)
}