# Nomenclatura IUPAC (PHP puro, sin base de datos) - imagen para Dokploy
FROM php:8.3-apache

# Módulos usados por los .htaccess del proyecto
RUN a2enmod headers expires rewrite \
 # Detrás de Traefik (Dokploy) el TLS termina en el proxy: marcar HTTPS
 # para que la cookie de sesión salga con el atributo Secure
 && printf 'SetEnvIf X-Forwarded-Proto "^https$" HTTPS=on\nServerName localhost\n' \
      > /etc/apache2/conf-available/proxy-https.conf \
 && a2enconf proxy-https \
 && mv "$PHP_INI_DIR/php.ini-production" "$PHP_INI_DIR/php.ini"

COPY --chown=www-data:www-data . /var/www/html/

EXPOSE 80
